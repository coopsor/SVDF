from __future__ import print_function

import pickle

from SVSignature import SignatureDeletion, SignatureInsertion, SignatureInversion, SignatureDuplicationTandem, SignatureTranslocation

def feature_record(alignment_current, alignment_next, ins_tra_flag=False, mode='general'):
    distance_on_read = alignment_next['q_start'] - alignment_current['q_end']
    if not alignment_current['is_reverse']:
        distance_on_reference = alignment_next['ref_start'] - alignment_current['ref_end']
        if alignment_next['is_reverse']:  # INV:+-
            if alignment_current['ref_end'] > alignment_next['ref_end']:
                distance_on_reference = alignment_next['ref_end'] - alignment_current['ref_start']
            else:
                distance_on_reference = alignment_current['ref_end'] - alignment_next['ref_start']
    else:
        distance_on_reference = alignment_current['ref_start'] - alignment_next['ref_end']
        if not alignment_next['is_reverse']:  # INV:-+
            if alignment_current['ref_end'] > alignment_next['ref_end']:
                distance_on_reference = alignment_next['ref_end'] - alignment_current['ref_start']
            else:
                distance_on_reference = alignment_current['ref_end'] - alignment_next['ref_start']
    deviation = distance_on_read - distance_on_reference
    chr_ = 1 if alignment_current['ref_id'] == alignment_next['ref_id'] else 0
    orientation = 1 if alignment_current['is_reverse'] == alignment_next['is_reverse'] else 0
    if mode == 'general':
        soft_clipped_length = 0
        for operation, length in alignment_current['cigar_tuple']:
            if operation == 4:
                soft_clipped_length += length
        soft_clipped_ratio = soft_clipped_length / alignment_current['infer_read_length']
        atgc_seq = alignment_current['atgc_seq']
        a_freq = atgc_seq.count('A') / len(atgc_seq)
        t_freq = atgc_seq.count('T') / len(atgc_seq)
        g_freq = atgc_seq.count('G') / len(atgc_seq)
        c_freq = atgc_seq.count('C') / len(atgc_seq)
        dup_freq = (atgc_seq.count('AAAAA') + atgc_seq.count('TTTTT') + atgc_seq.count('GGGGG') + atgc_seq.count(
            'CCCCC')) / len(atgc_seq)
        feature_data = {
            'chr': chr_,
            'orientation': orientation,
            'dis_read': distance_on_read,
            'dis_ref': distance_on_reference,
            'deviation': deviation,
            'map_quality': (alignment_current['mapping_quality'] + alignment_next['mapping_quality']) / 2,
            'soft_clipped_ratio': soft_clipped_ratio,
            'FLAG': alignment_current['FLAG'],
            'infer_read_length': alignment_current['infer_read_length'],
            'bin': alignment_current['bin'],
            'A': alignment_current['q_start'],
            'B': alignment_current['q_end'],
            'C': alignment_current['ref_id'],
            'D': alignment_current['ref_start'],
            'E': alignment_current['ref_end'],
            'soft_clipped_length': soft_clipped_length,
            'NM': alignment_current['NM'] if alignment_current['NM'] else 0,
            'a_freq': a_freq,
            't_freq': t_freq,
            'g_freq': g_freq,
            'c_freq': c_freq,
            'dup_freq': dup_freq,
        }
        return [feature_data, (alignment_current, alignment_next,  ins_tra_flag)]
    else:
        return (alignment_current, alignment_next, chr_, orientation, distance_on_read, distance_on_reference,
                deviation, ins_tra_flag)

def feature_read_segement(primary, supplementaries, mode = 'general'):
    read_name = primary.query_name
    alignments = [primary] + supplementaries
    alignment_list = []
    sg_list = []
    for alignment in alignments:
        if alignment.is_reverse:
            q_start = alignment.infer_read_length() - alignment.query_alignment_end
            q_end = alignment.infer_read_length() - alignment.query_alignment_start
        else:
            q_start = alignment.query_alignment_start
            q_end = alignment.query_alignment_end

        new_alignment_dict = {'read_name': read_name,
                              'q_start': q_start,
                              'q_end': q_end,
                              'ref_id': alignment.reference_id,
                              'ref_start': alignment.reference_start,
                              'ref_end': alignment.reference_end,
                              'is_reverse': alignment.is_reverse,
                              'mapping_quality': alignment.mapping_quality,
                              'cigar_tuple': alignment.cigartuples,
                              'FLAG': alignment.flag,
                              'infer_read_length': alignment.infer_read_length(),
                              'NM': alignment.get_tag('NM'),
                              'bin': alignment.bin,
                              'atgc_seq': alignment.query_sequence,
                              }
        alignment_list.append(new_alignment_dict)

    sorted_alignment_list = sorted(alignment_list, key=lambda aln: (aln['q_start'], aln['q_end']))
    for index in range(len(sorted_alignment_list) - 1):
        sg_list.append(feature_record(sorted_alignment_list[index], sorted_alignment_list[index + 1], mode=mode))
    if len(alignment_list) >= 3 and sorted_alignment_list[0]['ref_id'] != sorted_alignment_list[1]['ref_id']:
        sg_list.append(feature_record(sorted_alignment_list[0], sorted_alignment_list[-1], ins_tra_flag=True, mode=mode))

    return sg_list

def analyze_read_segments(bam, segement_data, options, mode = 'general'):
    sv_signatures = []
    if mode == 'general':
        feature_list = [list(sublist[0].values()) for sublist in segement_data]
        sig_list = [sublist[1] for sublist in segement_data]
        model = pickle.load(open('model/model_sim_indel_rf_test_.pkl', 'rb'))
        pred_type_list = model.predict(feature_list)

        for pred_type, feature_data, split_data in zip(pred_type_list, feature_list, sig_list):
            alignment_current = split_data[0]
            alignment_next = split_data[1]
            read_name = alignment_current['read_name']
            ref_chr = bam.getrname(alignment_current['ref_id'])
            orientation = feature_data[1]
            deviation = feature_data[4]
            if pred_type == 1 or split_data[2]:  # INS
                if not alignment_current['is_reverse']:
                    start = (alignment_current['ref_end'] + alignment_next['ref_start']) // 2 if not split_data[2] else min(alignment_current['ref_end'], alignment_next['ref_start'])
                else:
                    start = (alignment_current['ref_start'] + alignment_next['ref_end']) // 2 if not split_data[2] else min(alignment_current['ref_start'], alignment_next['ref_end'])
                end = start + deviation
                if end - start < options.min_sv_size:
                    continue
                sv_sig = (ref_chr, start, end, "suppl", read_name)
                sv_signatures.append(SignatureInsertion(*sv_sig))
            elif pred_type == 0:  # DEL
                if not alignment_current['is_reverse']:
                    start = alignment_current['ref_end']
                else:
                    start = alignment_next['ref_end']
                end = start - deviation
                if end - start < options.min_sv_size:
                    continue
                sv_sig = (ref_chr, start, end, "suppl", read_name)
                sv_signatures.append(SignatureDeletion(*sv_sig))
            elif pred_type == 2:
                # if distance_on_reference <= -options.min_sv_size:  # DUP: 重复区域的长度大于阈值且大于sv最小值
                if not alignment_current['is_reverse']:
                    start = alignment_next['ref_start']
                    end = alignment_current['ref_end']
                else:
                    start = alignment_current['ref_start']
                    end = alignment_next['ref_end']
                if end - start < options.min_sv_size:
                    continue
                sv_sig = (ref_chr, start, end, "suppl", read_name)
                sv_signatures.append(SignatureDuplicationTandem(*sv_sig))
            elif pred_type == 3:  # INV
                if not alignment_current['is_reverse']:  # +-
                    if alignment_next['ref_start'] - alignment_current['ref_end'] >= -options.segment_overlap_tolerance:
                        # if options.min_sv_size <= distance_on_reference <= options.max_sv_size:
                        start = alignment_current['ref_end']
                        end = alignment_next['ref_end']
                        sv_sig = (ref_chr, start, end, "suppl", read_name, "left_fwd")
                    elif alignment_current['ref_start'] - alignment_next['ref_end'] >= -options.segment_overlap_tolerance:
                        # if -options.max_sv_size <= distance_on_reference <= -options.min_sv_size:
                        start = alignment_next['ref_end']
                        end = alignment_current['ref_end']
                        sv_sig = (ref_chr, start, end, "suppl", read_name, "left_rev")
                    else:
                        continue
                else:  # -+
                    if alignment_next['ref_start'] - alignment_current['ref_end'] >= -options.segment_overlap_tolerance:
                        # if options.min_sv_size <= distance_on_reference <= options.max_sv_size:
                        start = alignment_current['ref_start']
                        end = alignment_next['ref_start']
                        sv_sig = (ref_chr, start, end, "suppl", read_name, "right_fwd")
                    elif alignment_current['ref_start'] - alignment_next['ref_end'] >= -options.segment_overlap_tolerance:
                        # if -options.max_sv_size <= distance_on_reference <= -options.min_sv_size:
                        start = alignment_next['ref_start']
                        end = alignment_current['ref_start']
                        sv_sig = (ref_chr, start, end, "suppl", read_name, "right_rev")
                    else:
                        continue
                if end - start < options.min_sv_size:
                    continue
                sv_signatures.append(SignatureInversion(*sv_sig))
            elif pred_type == 4:  # TRA
                ref_chr_next = bam.getrname(alignment_next['ref_id'])
                if orientation == 1:
                    # if distance_on_read <= options.segment_gap_tolerance:
                    if not alignment_current['is_reverse']:  # ++
                        if ref_chr < ref_chr_next:
                            start = alignment_current['ref_end']
                            end = alignment_next['ref_start']
                        else:
                            ref_chr, ref_chr_next = ref_chr_next, ref_chr
                            start = alignment_next['ref_start']
                            end = alignment_current['ref_end']
                        sv_sig = (ref_chr, start, 'fwd', ref_chr_next, end, 'fwd', "suppl",  read_name)
                    else:  # --
                        if ref_chr < ref_chr_next:
                            start = alignment_current['ref_start']
                            end = alignment_next['ref_end']
                        else:
                            ref_chr, ref_chr_next = ref_chr_next, ref_chr
                            start = alignment_next['ref_end']
                            end = alignment_current['ref_start']
                        sv_sig = (ref_chr, start, 'rev', ref_chr_next, end, 'rev', "suppl",  read_name)
                else:
                    # if distance_on_read <= options.segment_gap_tolerance:
                    if not alignment_current['is_reverse']:  # +-
                        if ref_chr < ref_chr_next:
                            start = alignment_current['ref_end']
                            end = alignment_next['ref_end']
                        else:
                            ref_chr, ref_chr_next = ref_chr_next, ref_chr
                            start = alignment_next['ref_end']
                            end = alignment_current['ref_end']
                        sv_sig = (ref_chr, start, 'fwd', ref_chr_next, end, 'rev', "suppl", read_name)
                    else:  # -+
                        if ref_chr < ref_chr_next:
                            start = alignment_current['ref_start']
                            end = alignment_next['ref_start']
                        else:
                            ref_chr, ref_chr_next = ref_chr_next, ref_chr
                            start = alignment_next['ref_start']
                            end = alignment_current['ref_start']
                        sv_sig = (ref_chr, start, 'rev', ref_chr_next, end, 'fwd', "suppl", read_name)

                sv_signatures.append(SignatureTranslocation(*sv_sig))

        return sv_signatures
    else:
        for sv_sig in segement_data:
            alignment_current = sv_sig[0]
            alignment_next = sv_sig[1]
            read_name = alignment_current['read_name']
            ref_chr = bam.getrname(alignment_current['ref_id'])
            chr_, orientation, distance_on_read, distance_on_reference, deviation, long_ins = sv_sig[2:]
            if chr_ == 1:
                if orientation == 1:
                    if distance_on_reference >= -options.min_sv_size or long_ins:  # INS
                        if deviation > 0:  # INS
                            if not alignment_current['is_reverse']:
                                start = (alignment_current['ref_end'] + alignment_next['ref_start']) // 2 if not long_ins else min(alignment_current['ref_end'], alignment_next['ref_start'])
                            else:
                                start = (alignment_current['ref_start'] + alignment_next['ref_end']) // 2 if not long_ins else min(alignment_current['ref_start'], alignment_next['ref_end'])
                            end = start + deviation
                            if end - start < options.min_sv_size:
                                continue
                            sv_sig = (ref_chr, start, end, "suppl", read_name)
                            sv_signatures.append(SignatureInsertion(*sv_sig))
                        elif deviation < 0:  # DEL
                            if not alignment_current['is_reverse']:
                                start = alignment_current['ref_end']
                            else:
                                start = alignment_next['ref_end']
                            end = start - deviation
                            sv_sig = (ref_chr, start, end, "suppl", read_name)
                            sv_signatures.append(SignatureDeletion(*sv_sig))
                        else:
                            continue
                    else:
                        # if distance_on_reference <= -options.min_sv_size:  #
                        if not alignment_current['is_reverse']:
                            start = alignment_next['ref_start']
                            end = alignment_current['ref_end']
                        else:
                            start = alignment_current['ref_start']
                            end = alignment_next['ref_end']
                        sv_sig = (ref_chr, start, end, "suppl", read_name)
                        sv_signatures.append(SignatureDuplicationTandem(*sv_sig))
                else:  # INV
                    if not alignment_current['is_reverse']:  # +-
                        if alignment_next['ref_start'] - alignment_current[
                            'ref_end'] >= -options.segment_overlap_tolerance:
                            # if options.min_sv_size <= distance_on_reference <= options.max_sv_size:
                            start = alignment_current['ref_end']
                            end = alignment_next['ref_end']
                            sv_sig = (ref_chr, start, end, "suppl", read_name, "left_fwd")
                        elif alignment_current['ref_start'] - alignment_next[
                            'ref_end'] >= -options.segment_overlap_tolerance:
                            # if -options.max_sv_size <= distance_on_reference <= -options.min_sv_size:
                            start = alignment_next['ref_end']
                            end = alignment_current['ref_end']
                            sv_sig = (ref_chr, start, end, "suppl", read_name, "left_rev")
                        else:
                            continue
                    else:  # -+
                        if alignment_next['ref_start'] - alignment_current[
                            'ref_end'] >= -options.segment_overlap_tolerance:
                            # if options.min_sv_size <= distance_on_reference <= options.max_sv_size:
                            start = alignment_current['ref_start']
                            end = alignment_next['ref_start']
                            sv_sig = (ref_chr, start, end, "suppl", read_name, "right_fwd")
                        elif alignment_current['ref_start'] - alignment_next[
                            'ref_end'] >= -options.segment_overlap_tolerance:
                            # if -options.max_sv_size <= distance_on_reference <= -options.min_sv_size:
                            start = alignment_next['ref_start']
                            end = alignment_current['ref_start']
                            sv_sig = (ref_chr, start, end, "suppl", read_name, "right_rev")
                        else:
                            continue
                    sv_signatures.append(SignatureInversion(*sv_sig))
            else:  # TRA
                ref_chr_next = bam.getrname(alignment_next['ref_id'])
                if orientation == 1:
                    # if distance_on_read <= options.segment_gap_tolerance:
                    if not alignment_current['is_reverse']:  # ++
                        if ref_chr < ref_chr_next:
                            start = alignment_current['ref_end']
                            end = alignment_next['ref_start']
                        else:
                            ref_chr, ref_chr_next = ref_chr_next, ref_chr
                            start = alignment_next['ref_start']
                            end = alignment_current['ref_end']
                        sv_sig = (ref_chr, start, 'fwd', ref_chr_next, end, 'fwd', "suppl", read_name)
                    else:  # --
                        if ref_chr < ref_chr_next:
                            start = alignment_current['ref_start']
                            end = alignment_next['ref_end']
                        else:
                            ref_chr, ref_chr_next = ref_chr_next, ref_chr
                            start = alignment_next['ref_end']
                            end = alignment_current['ref_start']
                        sv_sig = (ref_chr, start, 'rev', ref_chr_next, end, 'rev', "suppl", read_name)
                else:
                    # if distance_on_read <= options.segment_gap_tolerance:
                    if not alignment_current['is_reverse']:  # +-
                        if ref_chr < ref_chr_next:
                            start = alignment_current['ref_end']
                            end = alignment_next['ref_end']
                        else:
                            ref_chr, ref_chr_next = ref_chr_next, ref_chr
                            start = alignment_next['ref_end']
                            end = alignment_current['ref_end']
                        sv_sig = (ref_chr, start, 'fwd', ref_chr_next, end, 'rev', "suppl", read_name)
                    else:  # -+
                        if ref_chr < ref_chr_next:
                            start = alignment_current['ref_start']
                            end = alignment_next['ref_start']
                        else:
                            ref_chr, ref_chr_next = ref_chr_next, ref_chr
                            start = alignment_next['ref_start']
                            end = alignment_current['ref_start']
                        sv_sig = (ref_chr, start, 'rev', ref_chr_next, end, 'fwd', "suppl", read_name)
                sv_signatures.append(SignatureTranslocation(*sv_sig))

    return sv_signatures

