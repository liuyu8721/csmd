#!/bin/bash -l

# need to obtain from license
HUMAN_LCR=/public/users/liuyu/glue/db/LCR/RepBase.ref

pathoscope_python=/public/users/liuyu/csmd_test/code/csmd-master/pathoscope2/Clinical_PathoScope/pathoscope/pathoscope.py

module load bwa/v0.7.15
module load samtools/1.4
module load bedtools/v2.26.0
module load seqtk/master
module load RepeatMasker/4.0.7
module load bowtie2/2.3.3.1
module load R/3.3.3
module load python/2.7.13
module load blast/2.6.0+
export PATH=/public/software/fastq-tools-0.8/bin:$PATH
