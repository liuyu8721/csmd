# All rights reserved by Liuyu
# Author: Liuyu
#!/bin/sh
source ${CSMD_HOME}/configs.sh

PF_WORKDIR_P2_REFNAME=$1
PF_WORKDIR_P2_THREAD=$2
PF_WORKDIR_P2_INPUT=$3
PF_WORKDIR_P2_SAM=$4
PF_WORKDIR_P2_REPORT=$5

echo "#--------------------------------------------------------------#"
echo "Per-sample microbiome profiling"

if [ ! -f "${PF_WORKDIR_P2_INPUT}" ];then
  echo -e "\033[44;37;5m ERROR: \033[0m ${PF_WORKDIR_P2_INPUT} not exist"
  exit
fi

echo "Alignment to csmdSpecies database"
bowtie2 \
-p ${PF_WORKDIR_P2_THREAD} -k 10 --very-sensitive-local \
-x ${PF_WORKDIR_P2_REFNAME} \
-f ${PF_WORKDIR_P2_INPUT} \
-S csmd_temp.sam

python ${pathoscope_python} \
-t sam -e csmd -f csmd_temp.sam

rm csmd_temp.sam
mv updated_csmd_temp.sam ${PF_WORKDIR_P2_SAM}
mv csmd-sam-report.tsv ${PF_WORKDIR_P2_REPORT}

