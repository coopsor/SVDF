import os
import re


def main(bam_file, working_file):
    def get_f1(svtype, k):

        os.system("grep -E \'#|SVTYPE=INS|SVTYPE=DEL\' {0}/{0}.vcf > variants.vcf".format(working_file))
        os.system("/home/user/miniconda3/envs/tf-2/bin/bgzip variants.vcf")
        os.system("/home/user/miniconda3/envs/tf-2/bin/tabix variants.vcf.gz")
        result = subprocess.run(
            "/home/user/miniconda3/envs/tf-2/bin/truvari bench  -f /home/public_data/sv/long_read/CHM13/ref/hg38.fa -b /home/public_data/sv/long_read/CHM13/CHM13.vcf.gz -o svim_eval  -r {1} --passonly --sizemin 50 --sizefilt 50 -p 0.00 -c variants.vcf.gz ".format(
                svtype, k, ), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cmd_out = result.stderr.decode()

        try:
            precision = re.findall(r'"precision":\s+(\d+\.\d+)', cmd_out)[-1]
            recall = re.findall(r'"recall":\s+(\d+\.\d+)', cmd_out)[-1]
            f1 = re.findall(r'"f1":\s+(\d+\.\d+)', cmd_out)[-1]
            tp = re.findall(r'"TP-base":\s+(\d+)', cmd_out)[-1]
            final_out = precision + '\t' + recall + '\t' + f1 + '\n' + tp
            print(final_out)
            os.system("rm variants.vcf.gz")
            os.system("rm variants.vcf.gz.tbi")
            os.system("rm -r svim_eval")
        except:
            print(cmd_out)
            os.system("rm variants.vcf.gz")
            os.system("rm variants.vcf.gz.tbi")
            os.system("rm -r svim_eval")

        return final_out, tp

    # os.system(
    #     "/home/user/miniconda3/envs/tf-2/bin/python svdf.py alignment --read_type ONT  -s 3 {0}  {1} /home/public_data/sv/long_read/CHM13/ref/hg38.fa".format(
    #          working_file, bam_file))
    os.system("grep \'#\' {0}/variants.vcf > {0}/{0}.vcf ".format(working_file))
    os.system('grep -v \'#\' {0}/variants.vcf | grep -v \'0/0\' | grep -v \'\\\./\\\.\' >> {0}/{0}.vcf'.format(working_file))
    shift_list = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 500]
    for j in ['del', 'ins']:
        tp_list = []
        print(j)
        for k in shift_list:
            final_out, tp = get_f1(j, k)
            tp_list.append(tp)
        print(tp_list)
        sub_value = [int(tp_list[i + 1]) - int(tp_list[i]) for i in range(len(tp_list) - 1)]
        for i, value in enumerate(sub_value):
            print(value)

if __name__ == '__main__':
    test_files = [
        '/home/public_data/sv/long_read/CHM13/ONT/chm13_ont.bam',
    ]
    working_file_list = ['chm13_ont']
    for test_file, working_file in zip(test_files, working_file_list, ):
        main(test_file, working_file)
