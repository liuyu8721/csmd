# All rights reserved by Liuyu
# Author: Liuyu
#!/bin/sh
source ${CSMD_HOME}/configs.sh

CS_P4_SAMPLE=cs_p4_sample_$(date +%N)

CS_WORKDIR_P4_REFNAME=$1

CS_WORKDIR_P4_THREAD=$2

CS_P4_INPUT=$3
CS_P4_OUTPUT=$4

CS_P4_INPATH="$(dirname "${CS_P4_INPUT}")"
CS_P4_OUTPATH="$(dirname "${CS_P4_OUTPUT}")"

echo "#--------------------------------------------------------------#"
echo "CS Phase IV: Extra human genome (EHG) removal"

echo "CS Phase IVa: Align the remaining unmapped reads to EHG using Bowtie2 with --very-sensitive-local parameter"

mkdir -p ${CS_P4_OUTPATH}
CS_P4_RUNNING=${CS_P4_OUTPATH}/.cs_p4_tmp_$(date +%N)
rm -rf ${CS_P4_RUNNING}
mkdir -p ${CS_P4_RUNNING}

### Aligning the remaining reads to human sequence database with Bowtie2
bowtie2 -p ${CS_WORKDIR_P4_THREAD} --very-sensitive-local \
	-x ${CS_WORKDIR_P4_REFNAME} \
	-f ${CS_P4_INPUT} \
  -S ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.sam

samtools view -bS -f 4 ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.sam > ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.unmapped.bam
bamToFastq -i ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.unmapped.bam -fq ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.unmapped.fastq
module load seqtk/master
seqtk seq -A ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.unmapped.fastq > ${CS_P4_OUTPUT}

echo "CS Phase IVb: Removing intermediate files"
rm -rf ${CS_P4_RUNNING}
