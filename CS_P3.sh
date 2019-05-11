# All rights reserved by Liuyu
# Author: Liuyu
#!/bin/sh
source ${CSMD_HOME}/configs.sh

#CS_P3_SAMPLE=SimData10
CS_P3_SAMPLE=cs_p3_sample_$(date +%N)

CS_WORKDIR_P3_REFNAME=$1

CS_WORKDIR_P3_THREAD=$2

CS_P3_INPUT=$3

CS_P3_OUTPUT=$4

CS_P3_INPATH="$(dirname "${CS_P3_INPUT}")"
CS_P3_OUTPATH="$(dirname "${CS_P3_OUTPUT}")"

echo "#--------------------------------------------------------------#"
echo "CS Phase III: Low complexity reads (LCR) removal"

echo "CS Phase IIIa: Running RepeatMasker"

if [ -f "${CS_P3_OUTPUT}" ];then
  echo "Error: file ${CS_P3_OUTPUT} already existed"
  exit
else
  mkdir -p ${CS_P3_OUTPATH}
fi

CS_P3_RUNNING=${CS_P3_OUTPATH}/.cs_p3_tmp_$(date +%N)
rm -rf ${CS_P3_RUNNING}
mkdir -p ${CS_P3_RUNNING}

### Fastq2Fasta
#  [LQZ] to fix module loading
module load seqtk/master
seqtk seq -A ${CS_P3_INPUT} > ${CS_P3_RUNNING}/${CS_P3_SAMPLE}.hg38Removal.AHGRemoval.fasta

### Identify and then discard the low complexity sequences in the data
RepeatMasker -e ncbi -lib ${CS_WORKDIR_P3_REFNAME} -qq -pa ${CS_WORKDIR_P3_THREAD} \
${CS_P3_RUNNING}/${CS_P3_SAMPLE}.hg38Removal.AHGRemoval.fasta -dir ${CS_P3_RUNNING}

### Extract the result file after LCR identification
java -jar ${CSMD_HOME}/WGSparser.jar parseRM \
${CS_P3_RUNNING}/${CS_P3_SAMPLE}.hg38Removal.AHGRemoval.fasta \
${CS_P3_OUTPUT}

echo "CS Phase IIIb: Removing intermediate files"

rm -rf ${CS_P3_RUNNING}