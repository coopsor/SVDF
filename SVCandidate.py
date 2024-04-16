import re
import time
from collections import defaultdict
import numpy as np

class Candidate:

    def __init__(self, contig, start, end, sv_types, type,  members, support_fraction=".", genotype="./.", ref_reads=None, alt_reads=None):
        self.contig = contig
        self.start = start
        self.end = end
        self.sv_types = sv_types
        self.type = type
        self.members = members
        self.score = len(members)
        self.support_fraction = support_fraction
        self.genotype = genotype
        self.ref_reads = ref_reads
        self.alt_reads = alt_reads

    def get_source(self):
        return (self.contig, self.start, self.end)

    def get_vcf_entry(self):
        contig, start, end = self.get_source()
        filters = []
        if self.genotype == "0/0":
            filters.append("hom_ref")
        if self.ref_reads != None and self.alt_reads != None:
            dp_string = str(self.ref_reads + self.alt_reads)
        else:
            dp_string = "."
        info_template = "SVTYPE={0};END={1};SVLEN={2};SUPPORT={3}"
        if self.type == 'INV' or self.type == 'DUP':
            info_string = info_template.format(self.type, end, end - start, sum(2 if i == 'suppl' else 1 for i in self.sv_types))
        else:  # DEL, INS
            info_string = info_template.format(self.type, end, start - end if self.type == "DEL" else end - start, self.score)
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
            chrom=contig,
            pos=start,
            id="PLACEHOLDERFORID",
            ref="N",
            alt="<" + self.type + ">",
            qual=self.score,
            filter="PASS" if len(filters) == 0 else ";".join(filters),
            info=info_string,
            format="GT:DP:AD",
            samples="{gt}:{dp}:{ref},{alt}".format(gt=self.genotype, dp=dp_string, ref=self.ref_reads if self.ref_reads != None else ".",
                                                   alt=self.alt_reads if self.alt_reads != None else "."))

class CandidateBreakend(Candidate):
    def __init__(self, source_contig, source_start, source_direction, dest_contig, dest_start, dest_direction, members,
                 sv_types, support_fraction=".", genotype="./.", ref_reads=None, alt_reads=None):
        self.contig = source_contig
        #0-based source of the translocation (first base before the translocation)
        self.start = source_start
        self.direction = source_direction
        self.dest_contig = dest_contig
        #0-based destination of the translocation (first base after the translocation)
        self.dest_start = dest_start
        self.dest_direction = dest_direction
        self.members = members
        self.score = len(members)
        self.sv_types = sv_types
        self.type = "BND"
        self.support_fraction = support_fraction
        self.genotype = genotype
        self.ref_reads = ref_reads
        self.alt_reads = alt_reads

    def get_source(self):
        return (self.contig, self.start)

    def get_destination(self):
        return (self.dest_contig, self.dest_start)

    def get_vcf_entry(self):
        source_contig, source_start = self.get_source()
        dest_contig, dest_start = self.get_destination()

        if (self.direction == 'fwd') and (self.dest_direction == 'fwd'):
            alt_string = "N[{contig}:{start}[".format(contig=dest_contig, start=dest_start)
        elif (self.direction == 'fwd') and (self.dest_direction == 'rev'):
            alt_string = "N]{contig}:{start}]".format(contig=dest_contig, start=dest_start)
        elif (self.direction == 'rev') and (self.dest_direction == 'rev'):
            alt_string = "]{contig}:{start}]N".format(contig=dest_contig, start=dest_start)
        elif (self.direction == 'rev') and (self.dest_direction == 'fwd'):
            alt_string = "[{contig}:{start}[N".format(contig=dest_contig, start=dest_start)
        filters = []
        if self.genotype == "0/0":
            filters.append("hom_ref")
        info_template = "SVTYPE={0};SUPPORT={1}"
        info_string = info_template.format(self.type, sum(2 if i == 'suppl' else 1 for i in self.sv_types))
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
            chrom=source_contig,
            pos=source_start,
            id="PLACEHOLDERFORID",
            ref="N",
            alt=alt_string,
            qual=self.score,
            filter="PASS" if len(filters) == 0 else ";".join(filters),
            info=info_string,
            format="GT:DP:AD",
            samples="{gt}:{dp}:{ref},{alt}".format(gt=self.genotype, dp='.', ref=".", alt="."))

    def get_vcf_entry_reverse(self):
        source_contig, source_start = self.get_destination()
        dest_contig, dest_start = self.get_source()
        if (self.direction == 'rev') and (self.dest_direction == 'rev'):
            alt_string = "N[{contig}:{start}[".format(contig=dest_contig, start=dest_start)
        elif (self.direction == 'fwd') and (self.dest_direction == 'rev'):
            alt_string = "N]{contig}:{start}]".format(contig=dest_contig, start=dest_start)
        elif (self.direction == 'fwd') and (self.dest_direction == 'fwd'):
            alt_string = "]{contig}:{start}]N".format(contig=dest_contig, start=dest_start)
        elif (self.direction == 'rev') and (self.dest_direction == 'fwd'):
            alt_string = "[{contig}:{start}[N".format(contig=dest_contig, start=dest_start)
        filters = []
        if self.genotype == "0/0":
            filters.append("hom_ref")
        info_template="SVTYPE={0};SUPPORT={1}"
        info_string = info_template.format(self.type, sum(2 if member=='suppl' else 1 for member in self.sv_types))
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
                    chrom=source_contig,
                    pos=source_start,
                    id="PLACEHOLDERFORID",
                    ref="N",
                    alt=alt_string,
                    qual=self.score,
                    filter="PASS" if len(filters) == 0 else ";".join(filters),
                    info=info_string,
                    format="GT:DP:AD",
                    samples="{gt}:{dp}:{ref},{alt}".format(gt=self.genotype, dp=".", ref=".", alt="."))

def consolidate_clusters_unilocal(clusters):
    """Consolidate clusters to a list of (type, contig, mean start, mean end, cluster size, members) tuples."""
    consolidated_clusters = []
    for index, cluster in enumerate(clusters):
        sv_types = [member.signature for member in cluster]
        members = [member.read for member in cluster]
        if cluster[0].type == 'BND':
            start = np.median([member.get_source()[1] for member in cluster])
            dest_start = np.median([member.get_destination()[1] for member in cluster])
            source_direction = max([member.source_direction for member in cluster], key=[member.source_direction for member in cluster].count)
            dest_direction = max([member.dest_direction for member in cluster], key=[member.dest_direction for member in cluster].count)
            consolidated_clusters.append(
                    CandidateBreakend(cluster[0].get_source()[0], int(round(start)), source_direction, cluster[0].get_destination()[0], int(round(dest_start)), dest_direction, members, sv_types))
        else:
            start = np.median([member.get_source()[1] for member in cluster])
            end = np.median([member.get_source()[2] - member.get_source()[1] for member in cluster]) + start
            consolidated_clusters.append(
                    Candidate(cluster[0].get_source()[0], int(round(start)), int(round(end)), sv_types, cluster[0].type, members))

    return consolidated_clusters

def sorted_nicely(vcf_entries):
    """ Sort the given vcf entries (in the form ((contig, start, end), vcf_string, sv_type)) in the way that humans expect.
        e.g. chr10 comes after chr2
        Algorithm adapted from https://blog.codinghorror.com/sorting-for-humans-natural-sort-order/"""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    tuple_key = lambda entry: (alphanum_key(str(entry[0][0])), entry[0][1], entry[0][2])
    return sorted(vcf_entries, key=tuple_key)

def write_final_vcf(deletion_candidates,
                    novel_insertion_candidates,
                    duplication_candidates,
                    inversion_candidates,
                    translation_candidates,
                    contig_names,
                    contig_lengths,
                    types_to_output,
                    options):
    vcf_output = open(options.working_dir+'/variants.vcf', 'w')

    # Write header lines
    print("##fileformat=VCFv4.2", file=vcf_output)
    print("##fileDate={0}".format(time.strftime("%Y-%m-%d|%I:%M:%S%p|%Z|%z")), file=vcf_output)
    for contig_name, contig_length in zip(contig_names, contig_lengths):
        print("##contig=<ID={0},length={1}>".format(contig_name, contig_length), file=vcf_output)
    if "DEL" in types_to_output:
        print("##ALT=<ID=DEL,Description=\"Deletion\">", file=vcf_output)
    if "INS" in types_to_output:
        print("##ALT=<ID=INS,Description=\"Insertion\">", file=vcf_output)
    print("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">", file=vcf_output)
    print("##INFO=<ID=CUTPASTE,Number=0,Type=Flag,Description=\"Genomic origin of interspersed duplication seems to be deleted\">", file=vcf_output)
    print("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">", file=vcf_output)
    print("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">", file=vcf_output)
    print("##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Number of reads supporting this variant\">", file=vcf_output)
    print("##FILTER=<ID=hom_ref,Description=\"Genotype is homozygous reference\">", file=vcf_output)
    print("##FILTER=<ID=not_fully_covered,Description=\"Tandem duplication is not fully covered by a single read\">", file=vcf_output)
    print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", file=vcf_output)
    print("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">", file=vcf_output)
    print("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth for each allele\">", file=vcf_output)
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample", file=vcf_output)

    # Prepare VCF entries depending on command-line parameters
    vcf_entries = []
    if "DEL" in types_to_output:
        for candidate in deletion_candidates:
            vcf_entries.append((candidate.get_source(), candidate.get_vcf_entry(), "DEL"))
    if "INS" in types_to_output:
        for candidate in novel_insertion_candidates:
            vcf_entries.append((candidate.get_source(), candidate.get_vcf_entry(), "INS"))
    if "DUP" in types_to_output:
        for candidate in duplication_candidates:
            vcf_entries.append((candidate.get_source(), candidate.get_vcf_entry(), "DUP"))
    if "INV" in types_to_output:
        for candidate in inversion_candidates:
            vcf_entries.append((candidate.get_source(), candidate.get_vcf_entry(), "INV"))
    if "BND" in types_to_output:
        for candidate in translation_candidates:
            vcf_entries.append(((candidate.get_source()[0], candidate.get_source()[1], candidate.get_source()[1] + 1), candidate.get_vcf_entry(), "BND"))
            # vcf_entries.append(((candidate.get_destination()[0], candidate.get_destination()[1], candidate.get_destination()[1] + 1), candidate.get_vcf_entry_reverse(), "BND"))

    # Sort and write entries to VCF
    svtype_counter = defaultdict(int)
    for source, entry, svtype in sorted_nicely(vcf_entries):
        variant_id = "svdf.{svtype}.{number}".format(svtype=svtype, number=svtype_counter[svtype] + 1)
        entry_with_id = entry.replace("PLACEHOLDERFORID", variant_id, 1)
        svtype_counter[svtype] += 1
        print(entry_with_id, file=vcf_output)

    vcf_output.close()
