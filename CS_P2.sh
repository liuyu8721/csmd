# All rights reserved by Liuyu
# Author: Liuyu
#!/bin/sh
source ${CSMD_HOME}/configs.sh

CS_P2_SAMPLE=cs_p2_sample_$(date +%N)

CS_WORKDIR_P2_REFNAME=$1

CS_WORKDIR_P2_THREAD=$2

CS_P2_INPUT=$3

CS_P2_OUTPUT=$4

CS_P2_INPATH="$(dirname "${CS_P2_INPUT}")"
CS_P2_OUTPATH="$(dirname "${CS_P2_OUTPUT}")"

echo "#--------------------------------------------------------------#"
echo "CS Phase II: Additional human genome (AHG) removal"

echo "CS Phase IIa: Align the remaining unmapped reads to additional human genomes with BWA"

mkdir -p ${CS_P2_OUTPATH}

CS_P2_RUNNING=${CS_P2_OUTPATH}/.cs_p2_tmp_$(date +%N)
rm -rf ${CS_P2_RUNNING}
mkdir -p ${CS_P2_RUNNING}

### construct bwa sai file
bwa aln -t ${CS_WORKDIR_P2_THREAD} -f ${CS_P2_RUNNING}/${SAMPLE}.AHG.sai ${CS_WORKDIR_P2_REFNAME} ${CS_P2_INPUT}
	
### sai2sam
bwa samse -f ${CS_P2_RUNNING}/${SAMPLE}.AHG.sam ${CS_WORKDIR_P2_REFNAME} ${CS_P2_RUNNING}/${SAMPLE}.AHG.sai ${CS_P2_INPUT}

echo "CS Phase IIb: AHG mapped reads removal"

samtools view -bS -f 4 ${CS_P2_RUNNING}/${SAMPLE}.AHG.sam > ${CS_P2_RUNNING}/${SAMPLE}.AHG.Unmapped.bam
samtools sort -n ${CS_P2_RUNNING}/${SAMPLE}.AHG.Unmapped.bam -o ${CS_P2_RUNNING}/${SAMPLE}.AHG.Unmapped.Sorted.bam
bamToFastq -i ${CS_P2_RUNNING}/${SAMPLE}.AHG.Unmapped.Sorted.bam -fq ${CS_P2_OUTPUT}

echo "CS Phase IIc: removing intermediate files"

rm -rf ${CS_P2_RUNNING}

