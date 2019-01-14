# All rights reserved by Liuyu
# Author: Liuyu
#!/bin/sh
source ${CSMD_HOME}/configs.sh

CS_P1_SAMPLE=cs_p1_sample_$(date +%N)

CS_WORKDIR_P1_REFNAME=$1

CS_WORKDIR_P1_THREAD=$2

CS_WORKDIR_P1_R1=$3

CS_WORKDIR_P1_R2=$4

CS_P1_OUTPUT=$5
CS_P1_OUTPATH="$(dirname "${CS_P1_OUTPUT}")"

echo "#--------------------------------------------------------------#"
echo "CS Phase I: Human reference genome (hg19) removal"

echo "CS Phase Ia: Aligning reads to hg19 with BWA"

mkdir -p ${CS_P1_OUTPATH}

CS_P1_RUNNING=${CS_P1_OUTPATH}/.cs_p1_tmp_$(date +%N)
rm -rf ${CS_P1_RUNNING}
mkdir -p ${CS_P1_RUNNING}

### construct bwa sai file
bwa aln -t ${CS_WORKDIR_P1_THREAD} -f ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.R1.sai ${CS_WORKDIR_P1_REFNAME} ${CS_WORKDIR_P1_R1}
bwa aln -t ${CS_WORKDIR_P1_THREAD} -f ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.R2.sai ${CS_WORKDIR_P1_REFNAME} ${CS_WORKDIR_P1_R2}

### sai2sam
bwa sampe -f ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.sam ${CS_WORKDIR_P1_REFNAME} \
${CS_P1_RUNNING}/${CS_P1_SAMPLE}.R1.sai ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.R2.sai \
${CS_WORKDIR_P1_R1} ${CS_WORKDIR_P1_R2}

echo "CS Phase Ib: hg19 mapped reads removal"

### Removing mapping reads
samtools view -bS -f 4 ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.sam > ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.Unmapped.bam

### Extract Read 1 from unmapped bam files
samtools view -b -f 64 ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.Unmapped.bam > ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.Unmapped.R1.bam
samtools sort -n ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.Unmapped.R1.bam -o ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.Unmapped.R1.Sorted.bam
bamToFastq -i ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.Unmapped.R1.Sorted.bam -fq ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.Unmapped.R1.fastq
sed -i '1~4 s/$/A/g' ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.Unmapped.R1.fastq

### Extract Read 2 from unmapped bam files
samtools view -b -f 128 ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.Unmapped.bam > ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.Unmapped.R2.bam
samtools sort -n ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.Unmapped.R2.bam -o ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.Unmapped.R2.Sorted.bam
bamToFastq -i ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.Unmapped.R2.Sorted.bam -fq ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.Unmapped.R2.fastq
sed -i '1~4 s/$/B/g' ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.Unmapped.R2.fastq

### Combine Read 1 and Read 2 data
cat ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.Unmapped.R1.fastq ${CS_P1_RUNNING}/${CS_P1_SAMPLE}.Unmapped.R2.fastq \
	> ${CS_P1_OUTPUT}

echo "CS Phase Ic: removing intermediate files"
rm -rf ${CS_P1_RUNNING}