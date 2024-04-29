import sys
import os
import logging
import argparse

def parse_arguments(arguments = sys.argv[1:]):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""SVDF - Structural Variant Detection Framework""")

    subparsers = parser.add_subparsers(help='modes', dest='sub')

    parser_bam = subparsers.add_parser('call',
                                        help='Detect SVs from an existing alignment')
    parser_bam.add_argument('bam_file',
                             type=str,
                             help='Coordinate-sorted and indexed BAM file with aligned long reads')
    parser_bam.add_argument('--working_dir',
                             type=os.path.abspath,
                             help='Working and output directory. \
                                   Existing files in the directory are overwritten. \
                                   If the directory does not exist, it is created.')
    parser_bam.add_argument('-t', '--num_threads',
                            type=int,
                            default=16,
                            help='Number of threads to use')
    parser_bam.add_argument('--read_type',
                            type=str,
                            default='ONT',
                            help='Platform type for sequencing data, either "CLR", "CCS" or "ONT" (default: %(default)s)')
    parser_bam.add_argument('-d', '--depth',
                            type=int,
                            help='Sequencing depth of this dataset')
    parser_bam.add_argument('-s', '--min_support',
                            type=int,
                            help='Minimal number of supporting reads for one SV event')
    parser_bam.add_argument('--mode',
                            type=str,
                            default='general',
                            help='mode to use general for use filter automatically, or use sensitive for more sensitive detection with no filter')
    parser_bam.add_argument('--types',
                            type=str,
                            default="DEL,INS,DUP,INV,BND",
                            help='SV types to include in output VCF (default: %(default)s). \
                                  Give a comma-separated list of SV types. The possible SV types are: DEL (deletions), \
                                  INS (novel insertions), INV (inversions), DUP:TANDEM (tandem duplications), \
                                  DUP:INT (interspersed duplications), BND (breakends).')
    parser_bam.add_argument('--ref',
                            type=str,
                            nargs='*',
                            default=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X',
                                     'Y', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                                     'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
                                     'chr20', 'chr21', 'chr22', 'chrX', 'chrY'],
                            help='Reference genome file that the long reads were aligned to (FASTA)')
    parser_bam.add_argument('--skip_genotype',
                            action='store_true',
                            help='Skip genotyping and only call SVs')
    parser_bam.add_argument('--min_mapq',
                                      type=int,
                                      default=20,
                                      help='Minimum mapping quality of reads to consider (default: %(default)s). \
                                            Reads with a lower mapping quality are ignored.')
    parser_bam.add_argument('--min_sv_size',
                                      type=int,
                                      default=40,
                                      help='Minimum SV size to detect (default: %(default)s). \
                                            SVDF can potentially detect events of any size but is limited by the \
                                            signal-to-noise ratio in the input alignments. That means that more \
                                            accurate reads and alignments enable the detection of smaller events. \
                                            For current PacBio or Nanopore data, we would recommend a minimum size \
                                            of 40bp or larger.')
    parser_bam.add_argument('--max_sv_size',
                                      type=int,
                                      default=100000,
                                      help='Maximum SV size to detect (default: %(default)s). \
                                              This parameter is used to distinguish long deletions (and inversions) from \
                                              translocations which cannot be distinguished from the alignment alone. \
                                              Split read segments mapping far apart on the reference could either \
                                              indicate a very long deletion (inversion) or a translocation breakpoint. \
                                              SVDF calls a translocation breakpoint if the mapping distance is larger \
                                              than this parameter and a deletion (or inversion) if it is smaller or equal.')
    parser_bam.add_argument('--segment_overlap_tolerance',
                                      type=int,
                                      default=5,
                                      help='Maximum tolerated overlap between adjacent alignment segments (default: %(default)s). \
                                            This parameter applies to overlaps on the reference and the read. Example: \
                                            Deletions are detected from two subsequent segments of a split read that are mapped \
                                            far apart from each other on the reference. The segment overlap tolerance determines \
                                            the maximum tolerated length of an overlap between both segments on the read. If the \
                                            overlap between the two segments on the read is larger than this value, no deletion is called.')

    parser_bam = subparsers.add_parser('test',
                                       help='Detect SVs from an existing signatures')
    parser_bam.add_argument('bam_file',
                            type=str,
                            help='Coordinate-sorted and indexed BAM file with aligned long reads')
    parser_bam.add_argument('--working_dir',
                            type=os.path.abspath,
                            help='Working and output directory. \
                                           Existing files in the directory are overwritten. \
                                           If the directory does not exist, it is created.')
    parser_bam.add_argument('-t', '--num_threads',
                            type=int,
                            default=16,
                            help='Number of threads to use')
    parser_bam.add_argument('--read_type',
                            type=str,
                            default='ONT',
                            help='Platform type for sequencing data, either "CLR", "CCS" or "ONT" (default: %(default)s)')
    parser_bam.add_argument('-d', '--depth',
                            type=int,
                            help='Sequencing depth of this dataset')
    parser_bam.add_argument('-s', '--min_support',
                            type=int,
                            help='Minimal number of supporting reads for one SV event')
    parser_bam.add_argument('--mode',
                            type=str,
                            default='general',
                            help='mode to use general for use filter automatically, or use sensitive for more sensitive detection with no filter')
    parser_bam.add_argument('--types',
                            type=str,
                            default="DEL,INS,INV,DUP,BND",
                            help='SV types to include in output VCF (default: %(default)s). \
                                  Give a comma-separated list of SV types. The possible SV types are: DEL (deletions), \
                                  INS (novel insertions), INV (inversions), DUP:TANDEM (tandem duplications), \
                                  DUP:INT (interspersed duplications), BND (breakends).')
    parser_bam.add_argument('--skip_genotype',
                            action='store_true',
                            help='Skip genotyping and only call SVs')
    parser_bam.add_argument('--min_mapq',
                            type=int,
                            default=20,
                            help='Minimum mapping quality of reads to consider (default: %(default)s). \
                                            Reads with a lower mapping quality are ignored.')
    parser_bam.add_argument('--min_sv_size',
                            type=int,
                            default=40,
                            help='Minimum SV size to detect (default: %(default)s). \
                                            SVDF can potentially detect events of any size but is limited by the \
                                            signal-to-noise ratio in the input alignments. That means that more \
                                            accurate reads and alignments enable the detection of smaller events. \
                                            For current PacBio or Nanopore data, we would recommend a minimum size \
                                            of 40bp or larger.')
    parser_bam.add_argument('--max_sv_size',
                            type=int,
                            default=100000,
                            help='Maximum SV size to detect (default: %(default)s). \
                                              This parameter is used to distinguish long deletions (and inversions) from \
                                              translocations which cannot be distinguished from the alignment alone. \
                                              Split read segments mapping far apart on the reference could either \
                                              indicate a very long deletion (inversion) or a translocation breakpoint. \
                                              SVDF calls a translocation breakpoint if the mapping distance is larger \
                                              than this parameter and a deletion (or inversion) if it is smaller or equal.')
    parser_bam.add_argument('--segment_overlap_tolerance',
                            type=int,
                            default=5,
                            help='Maximum tolerated overlap between adjacent alignment segments (default: %(default)s). \
                                            This parameter applies to overlaps on the reference and the read. Example: \
                                            Deletions are detected from two subsequent segments of a split read that are mapped \
                                            far apart from each other on the reference. The segment overlap tolerance determines \
                                            the maximum tolerated length of an overlap between both segments on the read. If the \
                                            overlap between the two segments on the read is larger than this value, no deletion is called.')

    return parser.parse_args(arguments)
