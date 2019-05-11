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
cp ${CS_P4_INPUT} ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.input.fasta
for ref in `ls ${CS_WORKDIR_P4_REFNAME}.??.fasta.1.bt2*`
do
    m=$(echo ${ref} | awk -F "." '{print $(NF-3)}')
    echo ${m}
    bowtie2 -p ${CS_WORKDIR_P4_THREAD} --very-sensitive-local \
        -x ${CS_WORKDIR_P4_REFNAME}.${m}.fasta \
	-f ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.input.fasta \
        -S ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.${m}.sam

    rm ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.input.fasta

    samtools view -bS -f 4 ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.${m}.sam \
      > ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.${m}.unmapped.bam
    bamToFastq -i ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.${m}.unmapped.bam \
      -fq ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.${m}.unmapped.fastq
    seqtk seq -A ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.${m}.unmapped.fastq \
      > ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.${m}.fasta

    cp ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.${m}.fasta ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.input.fasta

    rm ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.${m}.sam
    rm ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.${m}.unmapped.bam
    rm ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.${m}.unmapped.fastq
done
cp ${CS_P4_RUNNING}/${CS_P4_SAMPLE}.input.fasta ${CS_P4_OUTPUT}

echo "CS Phase IVb: Removing intermediate files"
rm -rf ${CS_P4_RUNNING}
