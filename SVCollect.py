import itertools
import logging
import re
from collections import Counter

import pandas
import pysam
from SVIntra import analyze_alignment_indel
from SVInter import analyze_read_segments, feature_read_segement


def retrieve_other_alignments(main_alignment, bam):
    """Reconstruct other alignments of the same read for a given alignment from the SA tag"""
    # reconstructing other alignments from SA tag does not work if sequence of main_alignment is hard-clipped
    if main_alignment.get_cigar_stats()[0][5] > 0:
        return []
    try:
        sa_tag = main_alignment.get_tag("SA").split(";")
    except KeyError:
        return []
    other_alignments = []
    # For each other alignment encoded in the SA tag
    for element in sa_tag:
        # Read information from the tag
        fields = element.split(",")
        if len(fields) != 6:
            continue
        rname = fields[0]
        pos = int(fields[1])
        strand = fields[2]
        # CIGAR string encoded in SA tag is shortened
        cigar = fields[3]
        mapq = int(fields[4])
        nm = int(fields[5])

        # Generate an aligned segment from the information
        a = pysam.AlignedSegment()
        a.query_name = main_alignment.query_name
        a.query_sequence = main_alignment.query_sequence
        if strand == "+":
            a.flag = 2048
        else:
            a.flag = 2064
        a.reference_id = bam.get_tid(rname)
        a.reference_start = pos - 1
        try:
            a.mapping_quality = mapq
        except OverflowError:
            a.mapping_quality = 0
        a.cigarstring = cigar
        a.next_reference_id = -1
        a.next_reference_start = -1
        a.template_length = 0
        a.query_qualities = main_alignment.query_qualities
        a.set_tags([("NM", nm, "i")])

        other_alignments.append(a)
    return other_alignments

def read_vcf():
    # vcf_data = pandas.read_csv('/home/public_data/sv/simulate/pacbio_sim/sorted.vcf', sep='\t', comment='#', header=None)
    vcf_data = pandas.read_csv('/home/user/code/SVdetect/data/giab/HG002_SVs_Tier1_v0.6.vcf', sep='\t', comment='#', header=None)
    # vcf_data = pandas.read_csv('/home/public_data/sv/simulate/hack.vcf', sep='\t', comment='#', header=None)
    chr_values = vcf_data[0].tolist()

    start_values = vcf_data[1].tolist()
    svlen_values = []
    svtype_values = []
    for value, ori in zip(vcf_data[7], vcf_data[4]):
        try:
            svlen = abs(int(value.split("SVLEN=")[1].split(";")[0]))
            # svlen = abs(int(value.split("SVLEN=")[1]))
        except:
            svlen = None
        if value.split("SVTYPE=")[1].split(";")[0] == 'BND':
            matches = re.findall(r'\d+', ori)
            if len(matches) == 1:
                svlen = (str(22) if ori.split(':')[0][-1] == 'X' else str(23)) + ':' + matches[0] + ':' + str(
                    svlen)
            else:
                svlen = str(int(matches[0]) - 1) + ':' + matches[1] + ':' + str(svlen)
        svlen_values.append(svlen)
        svtype_values.append(value.split("SVTYPE=")[1].split(";")[0])
    truth_list = list(zip(start_values, svlen_values, svtype_values))
    ground_truth = []

    counter = Counter(chr_values)
    start_index = 0
    for element, count in counter.items():
        ground_truth.append(truth_list[start_index:start_index + count])
        start_index = start_index + count
    return ground_truth

def merge_cigar(sigs, data_intra, plat):
    i = 0
    max_merge = 500 if plat == 'CCS' else 500
    while i < len(sigs) - 1:
        diff = abs(sigs[i + 1].start - sigs[i].start) if sigs[i].type == 'INS' else abs(sigs[i + 1].start - sigs[i].end)
        if sigs[i].type == sigs[i + 1].type and diff <= max_merge:
            sigs[i].end = sigs[i].end + sigs[i + 1].end - sigs[i + 1].start
            sigs.pop(i + 1)
            data_intra.pop(i + 1)
        else:
            i += 1
    return sigs, data_intra

def analyze_alignment_file_coordsorted(contig, start, end, options):
    from svdf import bam
    sv_signatures, sv_signatures_inter = [], []
    data_list_intra = []
    segement_data = []
    for current_alignment in bam.fetch(contig, start, end):
        try:
            if current_alignment.is_unmapped or current_alignment.is_secondary or current_alignment.mapping_quality < options.min_mapq or current_alignment.reference_start<start:
                continue
            sigs, data_intra = analyze_alignment_indel(current_alignment, bam, current_alignment.query_name, options)

            if sigs:
                sigs, data_intra = merge_cigar(sigs, data_intra, options.read_type)
                sv_signatures.extend(sigs)
                data_list_intra.extend(data_intra)
            if not current_alignment.is_supplementary:
                supplementary_alignments = retrieve_other_alignments(current_alignment, bam)
                good_suppl_alns = [aln for aln in supplementary_alignments if
                                   not aln.is_unmapped and aln.mapping_quality >= options.min_mapq]
                sig_list = feature_read_segement(current_alignment, good_suppl_alns, options.mode)
                segement_data.extend(sig_list)

        except StopIteration:
            break
        except KeyboardInterrupt:
            logging.warning('Execution interrupted by user. Stop detection and continue with next step..')
            break
    for sig, data in zip(sv_signatures, data_list_intra):
        sig.data = data
    if segement_data:
        sv_signatures_inter = analyze_read_segments(bam, segement_data, options, options.mode)

    sv_signatures = sv_signatures + sv_signatures_inter
    return sv_signatures
