import os
import re


def calculate_metrics(input_string):
    pattern = r'Overall: (\d+) ([\d/ ]+)'
    match = re.search(pattern, input_string)

    if match:
        groups = match.groups()
        TP_FN_FP_temp = re.split(r'[ /]+', groups[1].strip())
        TP_FN_FP = TP_FN_FP_temp[:-1]

        # 将每五个数字作为一组，计算 Precision、Recall 和 F1 score
        metrics = []
        tp, fn, fp = 0, 0, 0
        for i in [0, 4]:
            TP = int(TP_FN_FP[i])
            FN = int(TP_FN_FP[i + 5])
            FP = int(TP_FN_FP[i + 10])

            # 计算 Precision、Recall 和 F1 score
            precision = TP / (TP + FP) if TP + FP > 0 else 0
            recall = TP / (TP + FN) if TP + FN > 0 else 0
            f1_score = 2 * (precision * recall) / (precision + recall) if precision + recall > 0 else 0
            tp, fn, fp = tp + TP, fn + FN, fp + FP
            print(str(precision) + '\t' + str(recall) + '\t' + str(f1_score))
        pre = tp / (tp + fp) if tp + fp > 0 else 0
        rec = tp / (tp + fn) if tp + fn > 0 else 0
        f1 = 2 * pre * rec / (pre + rec)
        print(str(pre) + '\t' + str(rec) + '\t' + str(f1))
        return metrics
    else:
        return None


def main(read_type, bam_file, working_file, support):
    # os.system(
    #     "/home/user/miniconda3/envs/tf-2/bin/python svdf.py test  --skip_genotype --read_type {0} -s {1}   {2} --min_sv_size 40  {3} ../data/hs37d5.fa".format(
    #         read_type, support, working_file, bam_file))
    # os.system("cat {0}/variants.vcf > {0}/{0}_{1}.vcf".format(working_file, support))

    file_list = ['tools/svdf', 'tools/svim', 'tools/sniffles', 'tools/cutesv']
    for i in file_list:
        print(i)
        result = subprocess.run(
            "/home/user/download/SURVIVOR-master/Debug/SURVIVOR eval  {0}/{0}_{1}.vcf /home/user/download/SURVIVOR-master/Debug/sim_12-24.bed 500 eval_res".format(
                working_file, support), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cmd_out = result.stdout.decode()
        calculate_metrics(cmd_out)


if __name__ == '__main__':
    test_files = [
        '/home/public_data/sv/sim_pbsim_v2/aln_clr_12-24.bam',
        '/home/public_data/sv/sim_pbsim_v2/aln_ont_12-24.bam',
    ]
    working_file_list = ['pbsim_clr_12-24', 'pbsim_ont_12-24']
    support_list = [3, 3]
    read_type_list = ['CLR', 'ONT']
    for read_type, test_file, working_file, support in zip(read_type_list, test_files, working_file_list, support_list):
        print(support)
        main(read_type, test_file, working_file, support)
