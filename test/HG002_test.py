import os
import re


def main(read_type, bam_file, working_file, support):
    os.system(
        "/home/user/miniconda3/envs/tf-2/bin/python svdf.py test --read_type {0}  -s {1}  {2}  {3} ../data/hs37d5.fa".format(
            read_type, support, working_file, bam_file))
    os.system("cat {0}/variants.vcf > {0}/{0}.vcf".format(working_file))
    os.system("grep \'#\' {0}/variants.vcf >{0}/{0}_gt.vcf".format(working_file))
    os.system("grep -v \'#\' {0}/variants.vcf | grep -v '0/0\|\./\.' >> {0}/{0}_gt.vcf".format(working_file))

    # os.system("svim alignment tools/svim {0} /data/huheng/dataset/human/HG002_PacBio_CCS_15kb/ref/hs37d5.fa --min_sv_size 30 ".format(bam_file))
    #     # os.system("cat tools/svim/variants.vcf >tools/svim/{0}_variants.vcf".format(file_type))
    #     # os.system("grep \'#\' tools/svim/variants.vcf >tools/svim/{0}_gt.vcf".format(file_type))
    #     # os.system("grep -v \'#\' tools/svim/variants.vcf | grep -v '0/0\|\./\.' >> tools/svim/{0}_gt.vcf".format(file_type))
    #
    #     # os.system("sniffles --input {0}  -t 16 --vcf tools/sniffles/variants.vcf --allow-overwrite ".format(bam_file, file_type))
    #     # os.system("cat tools/sniffles/variants.vcf > tools/sniffles/{0}.vcf ".format(file_type))
    #     # os.system("grep \'#\' tools/sniffles/variants.vcf > tools/sniffles/{0}_gt.vcf ".format(file_type))
    #     # os.system('grep -v \'#\' tools/sniffles/variants.vcf | grep -v \'0/0\' | grep -v \'\\\./\\\.\' >> tools/sniffles/{0}_gt.vcf'.format(
    #     #         file_type))
    #
    #     # if file_type in ['ont_12-24', 'ont_12-24_20x', 'ont_12-24_10x', 'ont_12-24_5x']:
    #     #             os.system("cuteSV {0} /data/huheng/dataset/human/HG002_PacBio_CCS_15kb/ref/hs37d5.fa tools/cutesv/variants.vcf ./ --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3  -s {1} --genotype".format(bam_file, support))
    #     # elif file_type in ['clr_12-24', 'clr_12-24_35x', 'clr_12-24_20x', 'clr_12-24_10x', 'clr_12-24_5x', 'HG003', 'HG004']:
    #     #     os.system("cuteSV {0} /data/huheng/dataset/human/HG002_PacBio_CCS_15kb/ref/hs37d5.fa tools/cutesv/variants.vcf ./ --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5 -s {1} --genotype".format(bam_file, support))
    #     # else:
    #     #     os.system("cuteSV {0} /data/huheng/dataset/human/HG002_PacBio_CCS_15kb/ref/hs37d5.fa tools/cutesv/variants.vcf ./ --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5  -s {1} --genotype".format(bam_file, support))
    #
    #     # os.system("cat tools/cutesv/variants.vcf > tools/cutesv/{0}.vcf ".format(file_type))
    #     # os.system("grep \'#\' tools/cutesv/variants.vcf > tools/cutesv/{0}_gt.vcf ".format(file_type))
    #     # os.system('grep -v \'#\' tools/cutesv/variants.vcf | grep -v \'0/0\' | grep -v \'\\\./\\\.\' >> tools/cutesv/{0}_gt.vcf'.format(file_type))
    #

    def get_f1(svtype):
        os.system(
            "grep -E \'#|SVTYPE=INS|SVTYPE=DEL\' {1}/{1}_gt.vcf > variants.vcf".format(svtype.upper(), working_file))
        # os.system("grep -E \'#|SVTYPE={0}\' {1}/{1}.vcf > variants.vcf".format(svtype.upper(), working_file))
        os.system("/home/user/miniconda3/envs/tf-2/bin/bgzip variants.vcf")
        os.system("/home/user/miniconda3/envs/tf-2/bin/tabix variants.vcf.gz")
        print('genotype')
        result = subprocess.run(
            "/home/user/miniconda3/envs/tf-2/bin/truvari bench  --gtcomp -f ./data/hs37d5.fa -b ./data/giab/12_22.vcf.gz  -o svim_eval --sizemin 50 --sizefilt 50 --passonly  -p 0.00 -c variants.vcf.gz   --includebed ./data/giab/HG002_SVs_Tier1_v0.6.bed".format(
                svtype), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cmd_out = result.stderr.decode()
        try:
            precision = re.findall(r'"precision":\s+(\d+\.\d+)', cmd_out)[-1]
            recall = re.findall(r'"recall":\s+(\d+\.\d+)', cmd_out)[-1]
            f1 = re.findall(r'"f1":\s+(\d+\.\d+)', cmd_out)[-1]
            final_out = precision + '\t' + recall + '\t' + f1 + '\n'
            print(final_out)
            os.system("rm variants.vcf.gz")
            os.system("rm variants.vcf.gz.tbi")
            os.system("rm -r svim_eval")
        except:
            print(cmd_out)
            os.system("rm variants.vcf.gz")
            os.system("rm variants.vcf.gz.tbi")
            os.system("rm -r svim_eval")

        return final_out

    with open('f1_genetype.txt', 'a+') as f:
        for j in ['indel', ]:
            final_out = get_f1(j)
            # line = j[1] + '\t' + working_file + '\t' + final_out
        #     f.write(line)
        # f.write('\n')
    print(bam_file + ' write successfully')


if __name__ == '__main__':
    test_files = [
        '/home/public_data/sv/long_read/clr_12-24/clr_12-24.bam',
        '/home/public_data/sv/long_read/clr_12-24/clr_12-24_35x.bam',
        '/home/public_data/sv/long_read/clr_12-24/clr_12-24_20x.bam',
        '/home/public_data/sv/long_read/clr_12-24/clr_12-24_10x.bam',
        '/home/public_data/sv/long_read/clr_12-24/clr_12-24_5x.bam',
        '/home/public_data/sv/long_read/ccs_12-24/ccs_12-24.bam',
        '/home/public_data/sv/long_read/ccs_12-24/ccs_12-24_10x.bam',
        '/home/public_data/sv/long_read/ccs_12-24/ccs_12-24_5x.bam',
        '/home/public_data/sv/long_read/ont_12-24/ont_12-24.bam',
        '/home/public_data/sv/long_read/ont_12-24/ont_12-24_20x.bam',
        '/home/public_data/sv/long_read/ont_12-24/ont_12-24_10x.bam',
        '/home/public_data/sv/long_read/ont_12-24/ont_12-24_5x.bam',
        '/home/public_data/sv/long_read/HG003/HG003_12_24.bam',
        '/home/public_data/sv/long_read/HG003/HG004_12_24.bam',
    ]
    working_file_list = [
        'clr', 'clr_35x', 'clr_20x', 'clr_10x', 'clr_5x', 'ccs', 'ccs_10x', 'ccs_5x', 'ont', 'ont_20x', 'ont_10x',
        'ont_5x', 'HG003', 'HG004'
    ]
    support_list = [10, 5, 4, 3, 2, 3, 2, 1, 10, 4, 3, 2, 5, 5]
    read_type_list = ['CLR', 'CLR', 'CLR', 'CLR', 'CLR', 'CCS', 'CCS', 'CCS', 'ONT', 'ONT', 'ONT', 'ONT', 'CLR', 'CLR']

    for read_type, test_file, working_file, support in zip(read_type_list, test_files, working_file_list, support_list):
        main(read_type, test_file, working_file, support)
