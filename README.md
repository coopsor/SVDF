# SVDF
SVDF is a long-read-based structural variant caller uses deep learning network.

SVDF is able to detect and genotype DEL/INS/DUP/INV/TRA with fast running speed.


Installation
------------

    #Install from github (requires Python 3.6.* or newer): installs all dependencies
    git clone https://github.com/hh22814/SVDF.git
    cd SVDF

Dependencies
------------
- *tensorflow=2.6.0* 
- *pandas*
- *numpy* 
- *pysam* 
- *numba*
- *scipy*

Input
-----

SVDF takes sorted and indexed alignment files in BAM format as inputs. And SVDF has been successfully tested on PacBio CLR, PacBio HiFi (CCS) and Oxford Nanopore data and alignment files produced by the read aligners `minimap2 <https://github.com/lh3/minimap2>`_, `pbmm2 <https://github.com/PacificBiosciences/pbmm2/>`_ , `NGMLR <https://github.com/philres/ngmlr>`.

Output
------

SVDF produces SV calls in the Variant Call Format (VCF).

Usage
----------------------
    python svdf.py call ./HG002_PB_70x_RG_HP10XtrioRTG.bam --working_dir ./ -s 10 -t 16  

    #'--working_dir' the work path of SVDF to store temporary data and output files, and the output files include called vcf file and the log file.
    #'-s' to filter the low-quality SVs. If you are not sure the parameter, you must specify the sequencing data type(--read_type) and read depth(--depth), SVDF will automatically calculate this value.
    #'-t' to specify the number of threads to use, the default is 16.
    #'--mode' to specify the mode of SVDF, the default is 'general', you can use 'sensitive' mode to call cancer genome SVs or achieve high recalls.
    #'--skip_genotype' to skip the genotype step, the default is False.
    #'--ref' to specify the chromosomes list to call SVs, the default is all chromosomes.

Contact
-------

If you experience any problems or have suggestions please create an issue or a pull request.

Citation
--------
Heng Hu, Runtian Gao, Wentao Gao, Bo Gao, Zhongjun Jiang, Murong Zhou, Guohua Wang, Tao Jiang, SVDF: enhancing structural variation detect from long-read sequencing via automatic filtering strategies, Briefings in Bioinformatics, Volume 25, Issue 4, July 2024, bbae336, https://doi.org/10.1093/bib/bbae336


License
-------

The project is licensed under the GNU General Public License.
