import os
import os
import subprocess
import re

import numpy as np
import pandas
import pandas as pd
import vcf


# def main(read_type, bam_file, working_file, support):
#     # os.system("pbsv discover {0} \'tools/pbsv/{1}.svsig.gz\' --tandem-repeats /data/huheng/dataset/human/HG002_PacBio_CCS_15kb/ref/hs37d5.trf.bed".format(bam_file, bam_file.split('/')[-1].split('.')[0]))
#     # os.system("pbsv call /data/huheng/dataset/human/HG002_PacBio_CCS_15kb/ref/hs37d5.fa tools/pbsv/{0}.svsig.gz tools/pbsv/{0}.vcf -t INS,DEL".format(bam_file.split('/')[-1].split('.')[0]))
#     #
#     os.system("/home/user/miniconda3/envs/tf-2/bin/python svdf.py test --read_type {0}  -s {1}  {2} --min_sv_size 40  {3} ../data/hs37d5.fa".format(read_type, support, working_file,bam_file))
#     os.system("cat {0}/variants.vcf > {0}/{0}.vcf".format(working_file))
#     os.system("grep \'#\' {0}/variants.vcf >{0}/{0}_gt.vcf".format(working_file))
#     os.system("grep -v \'#\' {0}/variants.vcf | grep -v '0/0\|\./\.' >> {0}/{0}_gt.vcf".format(working_file))
#
#     #
#     # # os.system("sniffles --input {0}  --vcf tools/sniffles/{1}.vcf --allow-overwrite".format(bam_file, bam_file.split('/')[-1].split('.')[0]))
#     # # os.system("cat <(cat tools/sniffles/{0}.vcf | grep \'^#\') <(cat tools/sniffles/{0}.vcf | grep -vE \'^#\' | grep \'DUP\|INS\|DEL\' | sed \'s/DUP/INS/g\' | sort -k1,1 -k2,2g)".format(bam_file.split('/')[-1].split('.')[0]))
#     # # if file_type in ['ont_12-24', 'ont_12-24_20x', 'ont_12-24_10x', 'ont_12-24_5x']:
#     # #     os.system("cuteSV {0} /data/huheng/dataset/human/HG002_PacBio_CCS_15kb/ref/hs37d5.fa tools/cutesv/variants.vcf ./ --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3  -s {1}".format(bam_file, support))
#     # # elif file_type in ['clr_12-24', 'clr_12-24_35x', 'clr_12-24_20x', 'clr_12-24_10x', 'clr_12-24_5x']:
#     # #     os.system("cuteSV {0} /data/huheng/dataset/human/HG002_PacBio_CCS_15kb/ref/hs37d5.fa tools/cutesv/variants.vcf ./ --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5 -s {1}".format(bam_file, support))
#     # # else:
#     # #     os.system("cuteSV {0} /data/huheng/dataset/human/HG002_PacBio_CCS_15kb/ref/hs37d5.fa tools/cutesv/variants.vcf ./ --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5  -s {1}".format(bam_file, support))
#     # # os.system("grep \'#\' tools/cutesv/variants.vcf > tools/cutesv/{0}.vcf ".format(bam_file.split('/')[-1].split('.')[0]))
#     # # os.system('grep -v \'#\' tools/cutesv/variants.vcf | grep -v \'0/0\' | grep -v \'\\\./\\\.\' >> tools/cutesv/{0}.vcf'.format(bam_file.split('/')[-1].split('.')[0]))
#     #
#     def get_f1(svtype):
#         # vcf_reader = vcf.Reader(open("{0}/{0}.vcf".format(working_file), 'r'))
#         # vcf_writer = vcf.Writer(open('variants_filtered.vcf', 'w'), vcf_reader)
#         # for record in vcf_reader:
#         #     rec_support = record.INFO.get('SUPPORT')
#         #     if int(rec_support) >= support:
#         #         vcf_writer.write_record(record)
#         # os.system("/where/to/install/bin/bcftools view -i \"FILTER==\'PASS\' && QUAL > {0}\" tools/svim/{1}.vcf > variants_filtered.vcf".format(support, bam_file.split('/')[-1].split('.')[0]))
#         os.system("grep -E \'#|SVTYPE=INS|SVTYPE=DEL\' {1}/{1}_gt.vcf > variants.vcf".format(svtype.upper(), working_file))
#         # os.system("grep -E \'#|SVTYPE={0}\' {1}/{1}.vcf > variants.vcf".format(svtype.upper(), working_file))
#         os.system("/home/user/miniconda3/envs/tf-2/bin/bgzip variants.vcf")
#         os.system("/home/user/miniconda3/envs/tf-2/bin/tabix variants.vcf.gz")
#         print('genotype')
#         result = subprocess.run(
#             "/home/user/miniconda3/envs/tf-2/bin/truvari bench  --gtcomp -f ./data/hs37d5.fa -b ./data/giab/12_22.vcf.gz  -o svim_eval --sizemin 50 --sizefilt 50 --passonly  -p 0.00 -c variants.vcf.gz   --includebed ./data/giab/HG002_SVs_Tier1_v0.6.bed".format(
#                 svtype), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         cmd_out = result.stderr.decode()
#         try:
#             precision = re.findall(r'"precision":\s+(\d+\.\d+)', cmd_out)[-1]
#             recall = re.findall(r'"recall":\s+(\d+\.\d+)', cmd_out)[-1]
#             f1 = re.findall(r'"f1":\s+(\d+\.\d+)', cmd_out)[-1]
#             final_out = precision + '\t' + recall + '\t' + f1 + '\n'
#             print(final_out)
#             os.system("rm variants.vcf.gz")
#             os.system("rm variants.vcf.gz.tbi")
#             os.system("rm -r svim_eval")
#         except:
#             print(cmd_out)
#             os.system("rm variants.vcf.gz")
#             os.system("rm variants.vcf.gz.tbi")
#             os.system("rm -r svim_eval")
#
#         return final_out
#
#     with open('f1_genetype.txt', 'a+') as f:
#         for j in ['del',]:
#             final_out = get_f1(j)
#             line = j[1] + '\t' + working_file + '\t' + final_out
#         #     f.write(line)
#         # f.write('\n')
#     print(bam_file + ' write successfully')
#
#
# if __name__ == '__main__':
#     test_files = [
#                 # '/home/public_data/sv/long_read/clr_12-24/clr_12-24.bam',
#                 #   '/home/public_data/sv/long_read/clr_12-24/clr_12-24_35x.bam',
#                 #   '/home/public_data/sv/long_read/clr_12-24/clr_12-24_20x.bam',
#                 #   '/home/public_data/sv/long_read/clr_12-24/clr_12-24_10x.bam',
#                 #   '/home/public_data/sv/long_read/clr_12-24/clr_12-24_5x.bam',
#                 #   '/home/public_data/sv/long_read/ccs_12-24/ccs_12-24.bam',
#                 #   '/home/public_data/sv/long_read/ccs_12-24/ccs_12-24_10x.bam',
#                 #   '/home/public_data/sv/long_read/ccs_12-24/ccs_12-24_5x.bam',
#                 #   '/home/public_data/sv/long_read/ont_12-24/ont_12-24.bam',
#                 #   '/home/public_data/sv/long_read/ont_12-24/ont_12-24_20x.bam',
#                 #   '/home/public_data/sv/long_read/ont_12-24/ont_12-24_10x.bam',
#                 #   '/home/public_data/sv/long_read/ont_12-24/ont_12-24_5x.bam',
#                 '/home/public_data/sv/long_read/HG003/HG003_12_24.bam',
#                 '/home/public_data/sv/long_read/HG003/HG004_12_24.bam',
#                   ]
#     working_file_list = [
#       # 'clr', 'clr_35x','clr_20x','clr_10x','clr_5x','ccs','ccs_10x','ccs_5x','ont', 'ont_20x', 'ont_10x', 'ont_5x',
#         'HG003', 'HG004'
#     ]
#     support_list = [5,5]
#     read_type_list = ['CLR','CLR']
#
#     for read_type, test_file, working_file, support in zip(read_type_list, test_files, working_file_list, support_list):
#         main(read_type, test_file, working_file, support)

# def main(bam_file, working_file, support):
#     # os.system("pbsv discover {0} \'tools/pbsv/{1}.svsig.gz\' --tandem-repeats /data/huheng/dataset/human/HG002_PacBio_CCS_15kb/ref/hs37d5.trf.bed".format(bam_file, bam_file.split('/')[-1].split('.')[0]))
#     # os.system("pbsv call /data/huheng/dataset/human/HG002_PacBio_CCS_15kb/ref/hs37d5.fa tools/pbsv/{0}.svsig.gz tools/pbsv/{0}.vcf -t INS,DEL".format(bam_file.split('/')[-1].split('.')[0]))
#
#     # os.system(
#     #     "/home/user/miniconda3/envs/tf-2/bin/python svdf.py alignment --read_type {0}  -s {1}  {2} --min_sv_size 40  {3} /home/public_data/sv/long_read/CHM13/ref/hg38.fa".format(
#     #         read, support, working_file, bam_file))
#     # os.system("cat {0}/variants.vcf | "
#     #           "sed 's/q6/PASS/g' > {0}/{0}.vcf".format(working_file))
#     # os.system("grep \'#\' {0}/variants.vcf > {0}/{0}.vcf ".format(working_file))
#     # os.system(
#     #     'grep -v \'#\' {0}/variants.vcf | grep -v \'0/0\' | grep -v \'\\\./\\\.\' >> {0}/{0}.vcf'.
#     #         format(working_file))
#
#     def get_f1(svtype,k):
#         # vcf_reader = vcf.Reader(open("{0}/{0}_gt.vcf".format(working_file), 'r'))
#         # vcf_writer = vcf.Writer(open('variants_filtered.vcf', 'w'), vcf_reader)
#         # for record in vcf_reader:
#         #     rec_support = record.INFO.get('SUPPORT')
#         #     if int(rec_support) >= support:
#         #         vcf_writer.write_record(record)
#         # os.system("/where/to/install/bin/bcftools view -i \"FILTER==\'PASS\' && QUAL > {0}\" tools/svim/{1}.vcf > variants_filtered.vcf".format(support, bam_file.split('/')[-1].split('.')[0]))
#         # os.system("cat variants_filtered.vcf > variants.vcf".format(svtype[0]))
#         # os.system("grep -E \'#|SVTYPE={0}\' {1}/{1}.vcf > variants.vcf".format(svtype.upper(), working_file))
#         os.system("grep -E \'#|SVTYPE=INS|SVTYPE=DEL\' {0}/{0}.vcf > variants.vcf".format(working_file))
#         os.system("/home/user/miniconda3/envs/tf-2/bin/bgzip variants.vcf")
#         os.system("/home/user/miniconda3/envs/tf-2/bin/tabix variants.vcf.gz")
#         result = subprocess.run(
#             "/home/user/miniconda3/envs/tf-2/bin/truvari bench  -f /home/public_data/sv/long_read/CHM13/ref/hg38.fa -b /home/public_data/sv/long_read/CHM13/CHM13.vcf.gz -o svim_eval_chm_  -r {1} --passonly --sizemin 50 --sizefilt 50 -p 0.00 -c variants.vcf.gz ".format(
#                 svtype,  k, ), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         cmd_out = result.stderr.decode()
#
#         try:
#             precision = re.findall(r'"precision":\s+(\d+\.\d+)', cmd_out)[-1]
#             recall = re.findall(r'"recall":\s+(\d+\.\d+)', cmd_out)[-1]
#             f1 = re.findall(r'"f1":\s+(\d+\.\d+)', cmd_out)[-1]
#             tp = re.findall(r'"TP-base":\s+(\d+)', cmd_out)[-1]
#             final_out = precision + '\t' + recall + '\t' + f1 + '\n'+tp
#             print(final_out)
#             os.system("rm variants.vcf.gz")
#             os.system("rm variants.vcf.gz.tbi")
#             # os.system("rm -r svim_eval_chm_")
#         except:
#             print(cmd_out)
#             os.system("rm variants.vcf.gz")
#             os.system("rm variants.vcf.gz.tbi")
#             os.system("rm -r svim_eval")
#
#         return final_out, tp
#
#
#     for j in ['del']:
#         tp_list = []
#         print(j)
#         for k in [0]:
#             final_out, tp = get_f1(j,k)
#             tp_list.append(tp)
#         print(tp_list)
#             # 计算相邻差值
#         sub_value = [int(tp_list[i + 1]) - int(tp_list[i]) for i in range(len(tp_list) - 1)]
#         for i, value in enumerate(sub_value):
#             print(value)
#         #     f.write(line)
#         # f.write('\n')
#     print(bam_file + ' write successfully')
#     return tp
#
# if __name__ == '__main__':
#     test_files = [
#
#         # '/home/public_data/sv/long_read/CHM13/CCS/chm13_ccs_v2.bam',
#         '/home/public_data/sv/long_read/CHM13/ONT/chm13_ont.bam',
#         '/home/public_data/sv/long_read/CHM13/ONT/chm13_ont.bam',
#         '/home/public_data/sv/long_read/CHM13/ONT/chm13_ont.bam',
#         '/home/public_data/sv/long_read/CHM13/ONT/chm13_ont.bam',
#         '/home/public_data/sv/long_read/CHM13/ONT/chm13_ont.bam',
#         '/home/public_data/sv/long_read/CHM13/ONT/chm13_ont.bam',
#         '/home/public_data/sv/long_read/CHM13/ONT/chm13_ont.bam',
#         '/home/public_data/sv/long_read/CHM13/ONT/chm13_ont.bam',
#         '/home/public_data/sv/long_read/CHM13/ONT/chm13_ont.bam',
#         '/home/public_data/sv/long_read/CHM13/ONT/chm13_ont.bam'
#     ]
#     working_file_list = ['chm13_ont','chm13_ont','chm13_ont','chm13_ont','chm13_ont','chm13_ont','chm13_ont','chm13_ont','chm13_ont','chm13_ont']
#     support_list = [500]
#     tp_list = []
#     for test_file, working_file, support in zip(test_files, working_file_list, support_list):
#         main(test_file, working_file, support)


def calculate_metrics(input_string):
    # 使用正则表达式提取字符串中的数据
    # print(input_string)
    # input_format = input_string.split('Missing SVs: ')[1]
    # input_format = input_format.split('\n')
    #
    # sort_input_format = sorted(input_format, key=lambda x: x.split('\t')[0])
    # for i in range(4, len(sort_input_format)):
    #     sort_input_format[i] = sort_input_format[i]+'\t'+str(int(sort_input_format[i].split(' ')[-1])-int(sort_input_format[i].split(' ')[2]))
    # print('\n'.join(sort_input_format))

    pattern = r'Overall: (\d+) ([\d/ ]+)'
    match = re.search(pattern, input_string)

    if match:
        groups = match.groups()
        TP_FN_FP_temp = re.split(r'[ /]+', groups[1].strip()) # 按空格分割字符串，得到每组数字
        #删除最后两个元素
        TP_FN_FP = TP_FN_FP_temp[:-1]

        # 将每五个数字作为一组，计算 Precision、Recall 和 F1 score
        metrics = []
        tp,fn,fp = 0, 0, 0
        for i in [0,4]:
            TP = int(TP_FN_FP[i])
            FN = int(TP_FN_FP[i + 5])
            FP = int(TP_FN_FP[i + 10])

            # 计算 Precision、Recall 和 F1 score
            precision = TP / (TP + FP) if TP + FP > 0 else 0
            recall = TP / (TP + FN) if TP + FN > 0 else 0
            f1_score = 2 * (precision * recall) / (precision + recall) if precision + recall > 0 else 0
            tp,fn,fp = tp+TP,fn+FN,fp+FP
            print(str(precision) + '\t' + str(recall) + '\t' + str(f1_score))
        pre = tp / (tp + fp) if tp + fp > 0 else 0
        rec = tp / (tp + fn) if tp+fn > 0 else 0
        f1 = 2 * pre * rec / (pre+rec)
        print(str(pre) + '\t' + str(rec) + '\t' + str(f1))
        return metrics
    else:
        return None


# def main(read_type, bam_file, working_file, support):
#
#     os.system("/home/user/miniconda3/envs/tf-2/bin/python svdf.py test  --skip_genotype --read_type {0} -s {1}   {2} --min_sv_size 40  {3} ../data/hs37d5.fa".format(read_type, support, working_file,bam_file))
#     os.system("cat {0}/variants.vcf | "
#     "sed 's/q5/PASS/g' > {0}/{0}_{1}.vcf".format(working_file,support))
#
#     result = subprocess.run(
#         "/home/user/download/SURVIVOR-master/Debug/SURVIVOR eval  {0}/{0}_{1}.vcf /home/user/download/SURVIVOR-master/Debug/sim_12-24.bed 500 eval_res".format(
#             working_file,support), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#     cmd_out = result.stdout.decode()
#
#     calculate_metrics(cmd_out)
#
#
# if __name__ == '__main__':
#     test_files = [
#         '/home/public_data/sv/sim_pbsim_v2/aln_clr_12-24.bam',
#         '/home/public_data/sv/sim_pbsim_v2/aln_ont_12-24.bam',
#         '/home/public_data/sv/sim_pbsim_v2/aln_ont_12-24.bam',
#         '/home/public_data/sv/sim_pbsim_v2/aln_ont_12-24.bam',
#         '/home/public_data/sv/sim_pbsim_v2/aln_ont_12-24.bam',
#         '/home/public_data/sv/sim_pbsim_v2/aln_ont_12-24.bam',
#         '/home/public_data/sv/sim_pbsim_v2/aln_ont_12-24.bam',
#         '/home/public_data/sv/sim_pbsim_v2/aln_ont_12-24.bam',
#         '/home/public_data/sv/sim_pbsim_v2/aln_ont_12-24.bam',
#                   ]
#     working_file_list = ['pbsim_clr_12-24','pbsim_clr_12-24','pbsim_clr_12-24','pbsim_clr_12-24','pbsim_clr_12-24','pbsim_clr_12-24','pbsim_clr_12-24',]
#     support_list = [3,4,5,6,7,8,9]
#     read_type_list = ['CLR','ONT','ONT','ONT','ONT','CLR','ONT','ONT','ONT','ONT',]
#     for read_type, test_file, working_file, support in zip(read_type_list, test_files, working_file_list, support_list):
#         print(support)
#         main(read_type, test_file, working_file, support)

