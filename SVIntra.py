from __future__ import print_function

from SVSignature import SignatureDeletion, SignatureInsertion


def analyze_cigar_indel(tuples, min_length):
    """Parses CIGAR tuples (op, len) and returns Indels with a length > minLength"""
    pos_ref = 0
    pos_read = 0
    indels = []
    soft_clipped_length = 0
    for operation, length in tuples:
        if operation == 0:  # alignment match
            pos_ref += length
            pos_read += length
        elif operation == 1:  # insertion
            if length >= min_length:
                indels.append((pos_ref, pos_read, length, "INS", operation))
            pos_read += length
        elif operation == 2:  # deletion
            if length >= min_length:
                indels.append((pos_ref, pos_read, length, "DEL", operation))
            pos_ref += length
        elif operation == 4:  # soft clip
            pos_read += length
            soft_clipped_length += length
        elif operation == 7 or operation == 8:  # match or mismatch
            pos_ref += length
            pos_read += length
    return indels, [soft_clipped_length, len(indels)]


def analyze_alignment_indel(alignment, bam, query_name, options):
    sv_signatures = []
    data_list = []
    ref_chr = bam.getrname(alignment.reference_id)
    ref_start = alignment.reference_start
    indels, cigar_fea = analyze_cigar_indel(alignment.cigartuples, options.min_sv_size)
    soft_clipped_length = cigar_fea[0]
    soft_clipped_ratio = soft_clipped_length / alignment.infer_read_length()
    indels_count = cigar_fea[1]
    total_mismatch = alignment.get_tag('NM') - indels_count if alignment.get_tag('NM') else 0
    atgc_seq = alignment.query_sequence
    a_freq = atgc_seq.count('A') / len(atgc_seq)
    t_freq = atgc_seq.count('T') / len(atgc_seq)
    g_freq = atgc_seq.count('G') / len(atgc_seq)
    c_freq = atgc_seq.count('C') / len(atgc_seq)
    dup_freq = (atgc_seq.count('AAAAA') + atgc_seq.count('TTTTT') + atgc_seq.count('GGGGG') + atgc_seq.count(
        'CCCCC')) / len(atgc_seq)

    for pos_ref, pos_read, length, typ, operation in indels:
        start = ref_start + pos_ref
        end = ref_start + pos_ref + length
        if typ == "DEL":
            sv_signatures.append(SignatureDeletion(ref_chr, start, end, "cigar", query_name))
        elif typ == "INS":
            sv_signatures.append(SignatureInsertion(ref_chr, start, end, "cigar", query_name))
        data = {'map_quality': int(alignment.mapping_quality),
                'soft_clipped_ratio': soft_clipped_ratio,
                'FLAG': alignment.flag,
                'infer_read_length': alignment.infer_read_length(),
                'indels_count': indels_count,  # alignment比对信息
                'bin': alignment.bin,  # cigar信息
                'start': start,
                'end': end,
                'pos_ref': pos_ref,
                'pos_read': pos_read,
                'length': length,
                'operation': operation,
                'A': alignment.query_alignment_start,
                'B': alignment.query_alignment_end,
                'C': alignment.reference_id,
                'D': alignment.reference_start,
                'E': alignment.reference_end,
                'soft': soft_clipped_length,
                # 'query_qualities': sum(alignment.query_qualities) / len(alignment.query_qualities),
                'NM': alignment.get_tag('NM'),
                'total_mismatch': total_mismatch,
                'a_freq': a_freq,  # 碱基信息
                't_freq': t_freq,
                'g_freq': g_freq,
                'c_freq': c_freq,
                'dup_freq': dup_freq,
                'label': 1
                }
        data_list.append(data)

    return sv_signatures, data_list
