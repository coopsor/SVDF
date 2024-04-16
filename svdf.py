import pandas
import sys
import os
import pickle
import logging
import time
from multiprocessing import Pool
import pysam
from time import strftime, localtime

from input_parsing import parse_arguments
from SVCollect import analyze_alignment_file_coordsorted
from SVCluster import form_bins, cluster_data
from SVFilter import filter_clusters
from SVCandidate import consolidate_clusters_unilocal, write_final_vcf
from SVGenotype import genotype

options = parse_arguments()
bam = pysam.AlignmentFile(options.bam_file, threads=options.num_threads)

def read_bed(truth_file):
    bed_data = pandas.read_csv(truth_file, sep='\t', header=None).values.tolist()
    bed_data = sorted(bed_data, key=lambda x: x[-1])
    chr_values, start_values, svlen_values, svtype_values = [], [], [], []
    for i in range(len(bed_data)):
        if bed_data[i][4] != 'TRA':
            chr_values.append(bed_data[i][0])
            start_values.append(bed_data[i][1])
            svlen_values.append(bed_data[i][3] - bed_data[i][1])
        else:
            chr_values.append(bed_data[i][0])
            start_values.append(bed_data[i][1])
            svlen_values.append(bed_data[i][2]+':'+str(bed_data[i][3]))
        svtype_values.append(bed_data[i][4])

    result_dict = {}
    for i in range(len(chr_values)):
        if chr_values[i] in result_dict:
            result_dict[chr_values[i]].append([start_values[i], svlen_values[i], svtype_values[i]])
        else:
            result_dict[chr_values[i]] = [[start_values[i], svlen_values[i], svtype_values[i]]]
    return result_dict

def multi_process(total_len, step, args=None):
    num_threads = options.num_threads
    analysis_pools = Pool(processes=int(num_threads))
    async_results = []
    chunk_size = total_len // num_threads
    for i in range(num_threads):
        start = i * chunk_size
        end = start + chunk_size if i < num_threads - 1 else total_len
        if step == 'collect':
            async_results.append(analysis_pools.starmap_async(analyze_alignment_file_coordsorted, [(args, start, end, options)]))
        elif step == 'cluster':
            async_results.append(analysis_pools.starmap_async(cluster_data, [(args[0][start:end], args[1])]))
        else:
            async_results.append(analysis_pools.starmap_async(genotype, [(args[0][start:end], args[1], options)]))

    analysis_pools.close()
    analysis_pools.join()
    results = []
    for async_result in async_results:
        result = async_result.get()
        results.extend(result)
    return [item for sublist in results for item in sublist]

def main():
    # Fetch command-line options
    if not options.sub:
        print("Please choose one of the two modes ('alignment' or 'test'). See --help for more information.")
        return

    # Set up logging
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-7.7s]  %(message)s")
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)

    # Create working dir if it does not exist
    if not os.path.exists(options.working_dir):
        os.makedirs(options.working_dir)

    # Create log file
    fileHandler = logging.FileHandler(
        "{0}/SVDF_{1}.log".format(options.working_dir, strftime("%y%m%d_%H%M%S", localtime())), mode="w")
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    logging.info("****************** Start SVDF ******************")
    logging.info("CMD: python3 {0}".format(" ".join(sys.argv)))
    logging.info("WORKING DIR: {0}".format(os.path.abspath(options.working_dir)))
    for arg in vars(options):
        logging.info("PARAMETER: {0}, VALUE: {1}".format(arg, getattr(options, arg)))

    if not options.min_support:
        try:
            options.min_support = options.depth // 10
            if options.read_type == 'ONT' or options.read_type == 'CLR':
                options.min_support = options.min_support + 2
        except:
            print(
                "Please specify the '-s' to filter the low-quality SVs. If you are not sure the parameter, you must specify the sequencing data type(--read_type) and depth(--depth).")
            return

    logging.info("****************** STEP 1: COLLECT ******************")
    if options.sub == 'call':
        logging.info("MODE: call")
        logging.info("INPUT: {0}".format(os.path.abspath(options.bam_file)))
        try:
            bam.check_index()
        except ValueError:
            logging.warning(
                "Input BAM file is missing a valid index. Please generate with 'samtools faidx'. Continuing without genotyping for now..")
        except AttributeError:
            logging.warning(
                "pysam's .check_index raised an Attribute error. Something is wrong with the input BAM file.")
            return
        sv_signatures = []
        ref_list = bam.get_index_statistics()
        for ref in ref_list:
            if ref.mapped == 0:
                continue
            if str(ref[0]) not in options.ref:
                continue
            ref_len = bam.get_reference_length(ref[0])
            logging.info("Processing ref {0}...".format(ref[0]))
            sv_signatures.extend(multi_process(ref_len, 'collect', ref[0]))
            logging.info("Processed ref {0}...".format(ref[0]))
        with open(options.working_dir + '/signatures_test.pkl', 'wb') as temp:
            pickle.dump(sv_signatures, temp)
    else:
        with open(options.working_dir+'/signatures_test.pkl', 'rb') as f:
            sv_signatures = pickle.load(f)

    deletion_signatures = [ev for ev in sv_signatures if ev.type == "DEL"]
    insertion_signatures = [ev for ev in sv_signatures if ev.type == "INS"]
    duplication_signatures = [ev for ev in sv_signatures if ev.type == "DUP"]
    inversion_signatures = [ev for ev in sv_signatures if ev.type == "INV"]
    translation_signatures = [ev for ev in sv_signatures if ev.type == "BND"]

    logging.info("Found {0} signatures for deleted regions.".format(len(deletion_signatures)))
    logging.info("Found {0} signatures for inserted regions.".format(len(insertion_signatures)))
    logging.info("Found {0} signatures for duplicated regions.".format(len(duplication_signatures)))
    logging.info("Found {0} signatures for inverted regions.".format(len(inversion_signatures)))
    logging.info("Found {0} signatures for translocated regions.".format(len(translation_signatures)))

    logging.info("****************** STEP 2: CLUSTER ******************")
    start_time = time.time()
    signature_clusters = []
    bin_depth_list = []
    for element_signature in [deletion_signatures, insertion_signatures, inversion_signatures, duplication_signatures, translation_signatures]:
        if not element_signature:
            continue
        signature_bin, bin_depth = form_bins(element_signature, 1000)
        bin_depth_list.append(bin_depth)
        signature_clusters.extend(multi_process(len(signature_bin), 'cluster', (signature_bin, bin_depth)))

    logging.info("Finished clustering..spend time: {0}".format(time.time() - start_time))

    logging.info("****************** STEP 3: FILTER ******************")
    start_time = time.time()
    del_min_support = options.min_support
    if options.read_type == 'ONT' or options.read_type == 'CLR':
        del_min_support = options.min_support - 1 if sum(bin_depth_list[0:2]) / 2 < 5 and bin_depth_list[0] > bin_depth_list[1] else options.min_support
        del_min_support = del_min_support if del_min_support > 1 else 2

    intra_cluster, inter_cluster = [], []
    for cluster in signature_clusters:
        sv_types = [obj.signature == 'cigar' for obj in cluster]
        if all(sv_types) and len(cluster) < 20:
            intra_cluster.append(cluster)
        else:
            factor = 2 if options.min_support >= 10 or (cluster[0].type != 'INS' and cluster[0].type != 'DEL') else 1
            filter_one = sum(factor if i is False else 1 for i in sv_types)
            if (cluster[0].type == 'DEL' and len(sv_types) >= del_min_support) or (cluster[0].type != 'DEL' and filter_one >= options.min_support):
                inter_cluster.append(cluster)
    if options.mode == 'general':
        intra_cluster = filter_clusters(intra_cluster)
    signature_clusters = []
    for cluster in intra_cluster:
        if (cluster[0].type == 'DEL' and len(cluster) >= del_min_support) or (cluster[0].type != 'DEL' and len(cluster) >= options.min_support):
            signature_clusters.append(cluster)
    signature_clusters.extend(inter_cluster)

    logging.info("Finished filtering.. spend time: {0}".format(time.time() - start_time))

    logging.info("****************** STEP 4: SVCALL ******************")
    signature_candidate = sorted(consolidate_clusters_unilocal(signature_clusters),
                                 key=lambda cluster: (cluster.contig, cluster.start))
    deletion_candidates = [i for i in signature_candidate if i.type == 'DEL']
    insertion_candidates = [i for i in signature_candidate if i.type == 'INS']
    duplication_candidates = [i for i in signature_candidate if i.type == 'DUP']
    inversion_candidates = [i for i in signature_candidate if i.type == 'INV']
    translation_candidates = [i for i in signature_candidate if i.type == 'BND']

    logging.info("Final deletion candidates: {0}".format(len(deletion_candidates)))
    logging.info("Final novel insertion candidates: {0}".format(len(insertion_candidates)))
    logging.info("Final duplication candidates: {0}".format(len(duplication_candidates)))
    logging.info("Final inversion candidates: {0}".format(len(inversion_candidates)))
    logging.info("Final translation candidates: {0}".format(len(translation_candidates)))

    # with open(options.working_dir+'/candidates_test.pkl', 'wb') as f:
    #     pickle.dump([deletion_candidates, insertion_candidates, duplication_candidates, inversion_candidates, translation_candidates], f)
    # with open(options.working_dir+'/candidates_test.pkl', 'rb') as f:
    #     cand = pickle.load(f)
    #     deletion_candidates, insertion_candidates, duplication_candidates, inversion_candidates, translation_candidates = cand
    if not options.skip_genotype:
        logging.info("****************** STEP 5: GENOTYPE ******************")
        logging.info("Genotyping deletions..")
        deletion_candidates = multi_process(len(deletion_candidates), 'genotype', (deletion_candidates, "DEL"))
        logging.info("Genotyping insertions..")
        insertion_candidates = multi_process(len(insertion_candidates), 'genotype', (insertion_candidates, "INS"))
        logging.info("Genotyping duplication..")
        duplication_candidates = multi_process(len(duplication_candidates), 'genotype', (duplication_candidates, "DUP"))
        logging.info("Genotyping inversions..")
        inversion_candidates = multi_process(len(inversion_candidates), 'genotype', (inversion_candidates, "INV"))
        logging.info("Genotyping translations..")
        translation_candidates = multi_process(len(translation_candidates), 'genotype', (translation_candidates, "BND"))

    types_to_output = [entry.strip() for entry in options.types.split(",")]
    write_final_vcf(deletion_candidates,
                    insertion_candidates,
                    duplication_candidates,
                    inversion_candidates,
                    translation_candidates,
                    bam.references,
                    bam.lengths,
                    types_to_output,
                    options)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        logging.error(e, exc_info=True)