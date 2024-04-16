import re
import pandas as pd

def convert_bed(excel_file):
    df = pd.read_excel(excel_file)
    #提取SV_id，chrom1,pos1,chrom2,pos2,svtype,SV_size or breakpoints distance,SV_type列
    df = df[['Chrom1', 'Pos1', 'Chrom2', 'Pos2', 'SV_size or breakpoints distance', 'SV_type']]
    df.loc[df['SV_type'] == 'DEL', 'SV_size or breakpoints distance'] = -df['SV_size or breakpoints distance']
    df.to_csv('tumor_v.bed', sep='\t', header=False, index=False)

# convert_bed('/home/user/download/SURVIVOR-master/Debug/13059_2022_2816_MOESM5_ESM.xlsx')

def read_bed(truth_file):
    bed_data = pd.read_csv(truth_file, sep='\t', header=None).values.tolist()
    # bed_data = sorted(bed_data, key=lambda x: x[-1])
    chr_values, start_values, svlen_values, svtype_values, sv_id = [], [], [], [], []
    for i in range(len(bed_data)):
        if bed_data[i][5] != 'TRA':
            chr_values.append(bed_data[i][0])
            start_values.append(bed_data[i][1])
            svlen_values.append(abs(bed_data[i][4]))
        else:
            chr_values.append(bed_data[i][0])
            start_values.append(bed_data[i][1])
            svlen_values.append(bed_data[i][2] + ':' + str(bed_data[i][3]))
        svtype_values.append(bed_data[i][5])
        sv_id.append(i)

    result_dict = {}
    for i in range(len(chr_values)):
        if chr_values[i] in result_dict:
            result_dict[chr_values[i]].append([start_values[i], svlen_values[i], svtype_values[i], sv_id[i]])
        else:
            result_dict[chr_values[i]] = [[start_values[i], svlen_values[i], svtype_values[i], sv_id[i]]]
    return result_dict


def main(read_type, bam_file, working_file, support):
    # os.system("/home/user/miniconda3/envs/tf-2/bin/python svdf.py test --mode sensitive  --read_type {0} -s {1} --skip_genotype  {2} --min_sv_size 50  {3} /home/public_data/sv/long_read/CHM13/ref/hg38.fa".format(read_type, support, working_file,bam_file))
    # os.system("cat {0}/variants.vcf > {0}/{0}.vcf".format(working_file))

    ground_truth = read_bed('/home/user/download/SURVIVOR-master/Debug/tumor_h.bed')
    result = {}
    with open("../{0}/{0}.vcf".format(working_file), 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            rec = line.split('\t')

            if rec[0] not in ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                              'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
                              'chr20', 'chr21', 'chr22']:
                continue
            chr, svstart = rec[0], int(rec[1])
            svtype = re.search(r'SVTYPE=([^;]+)', rec[7]).group(1)
            if svtype != 'BND':
                matches = re.search(r'SVLEN=([^;]+)', rec[7])
                if matches:
                    svlen = int(matches.group(1))
                else:
                    svlen = int(re.search(r'END=([^;]+)', rec[7]).group(1)) - svstart
                for truth in ground_truth[chr]:
                    truth_start, truth_len, truth_type, truth_id = truth[0], truth[1], truth[2], truth[3]
                    if svtype == truth_type:
                        if abs(truth_start - svstart) <= 1000 and 0.5 <= min(truth_len, abs(svlen)) / max(
                                truth_len, abs(svlen)):
                            if svtype in result:
                                result[svtype].append(truth_id)
                            else:
                                result[svtype] = [truth_id]
                            break
            else:
                matches = re.findall(r'chr(\d+):(\d+)', rec[4])
                if matches:
                    ori_chr, ori_pos = matches[0]
                    ori_chr = 'chr' + ori_chr
                else:
                    continue
                for truth in ground_truth[chr]:
                    truth_start, truth_end, truth_type, truth_id = truth[0], truth[1], truth[2], truth[3]
                    if 'TRA' == truth_type:
                        end_contig, end_pos = truth_end.split(':')[0], truth_end.split(':')[1]
                        if abs(truth_start - svstart) <= 1000 and abs(
                                float(end_pos) - float(ori_pos)) <= 1000 and end_contig == ori_chr:
                            if svtype in result:
                                result[svtype].append(truth_id)
                            else:
                                result[svtype] = [truth_id]
                            break
    return result

if __name__ == '__main__':
    test_files = [
        '/home/public_data/sv/long_read/tumor_clr.PBMM2.bam',
        '/home/public_data/sv/long_read/tumor.ont.bam',
    ]
    working_file_list = ['tumor_clr', 'tumor_ont']
    # working_file_list = ['pbsim_clr_12-24', 'pbsim_ont_12-24']
    support_list = [5, 3]
    read_type_list = ['CLR', 'ONT']
    result = []
    for read_type, test_file, working_file, support in zip(read_type_list, test_files, working_file_list, support_list):
        result.append(main(read_type, test_file, working_file, support))

    del_num = list(set(result[0]['DEL'] + result[1]['DEL']))
    ins_num = list(set(result[0]['INS'] + result[1]['INS']))
    dup_num = list(set(result[0]['DUP'] + result[1]['DUP']))
    ins_dup = ins_num + dup_num
    inv_num = list(set(result[0]['INV'] + result[1]['INV']))
    bnd_num = list(set(result[0]['BND'] + result[1]['BND']))
    max_length = max(len(del_num), len(ins_dup), len(inv_num), len(bnd_num),
                     sum([len(del_num), len(ins_dup), len(inv_num), len(bnd_num)]))

    df = pd.DataFrame({
        'del_num': del_num + [None] * (max_length - len(del_num)),
        'ins_dup': ins_dup + [None] * (max_length - len(ins_dup)),
        'inv_num': inv_num + [None] * (max_length - len(inv_num)),
        'bnd_num': bnd_num + [None] * (max_length - len(bnd_num)),
        'total': del_num + ins_dup + inv_num + bnd_num
    })

    df.fillna('', inplace=True)

    # 写入CSV文件
    df.to_csv('output_h.csv', index=False)
