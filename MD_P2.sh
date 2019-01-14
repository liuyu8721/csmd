# All rights reserved by Liuyu
# Author: Liuyu
#!/bin/sh
source ${CSMD_HOME}/configs.sh

MD_WORKDIR_P2_REFNAME=$1
MD_WORKDIR_P2_THREAD=$2
MD_WORKDIR_P2_INPUT=$3
MD_WORKDIR_P2_OUTDIR=$4

echo "#--------------------------------------------------------------#"
echo "MD Phase II: Species discovery"
if [ ! -f "${MD_WORKDIR_P2_INPUT}" ];then
  echo -e "\033[44;37;5m ERROR: \033[0m ${MD_WORKDIR_P2_INPUT} not exist"
  exit
fi

if [ ! -d "${MD_WORKDIR_P2_OUTDIR}" ];then
  echo -e "\033[44;37;5m ERROR:\033[0m OUTPATH can not be found"
  exit
fi

echo "MD Phase IIa: Aligned to SpeciesRepresentative DB"
> ${MD_WORKDIR_P2_OUTDIR}/Pooled_Bac_headers.sam
> ${MD_WORKDIR_P2_OUTDIR}/Pooled_Bac_noHeaders.sam

MD_WORKDIR_P2_REFLIST=(${MD_WORKDIR_P2_REFNAME}*.rev.1.bt2*)
for ref in "${MD_WORKDIR_P2_REFLIST[@]}"
do
  ref_prefix=$(echo ${ref} | awk -F '[.]' '{for(i=1;i<NF-2;i++){print $i}}')
  echo ${ref_prefix}
  bowtie2 \
    -p ${MD_WORKDIR_P2_THREAD} --very-sensitive -k 10 \
    -x ${ref_prefix} \
    -f ${MD_WORKDIR_P2_INPUT} \
    -S ${MD_WORKDIR_P2_OUTDIR}/Pooled_bac.sam

  samtools view -S -F 4 ${MD_WORKDIR_P2_OUTDIR}/Pooled_bac.sam >> ${MD_WORKDIR_P2_OUTDIR}/Pooled_Bac_noHeaders.sam
  samtools view -S -H ${MD_WORKDIR_P2_OUTDIR}/Pooled_bac.sam >> ${MD_WORKDIR_P2_OUTDIR}/Pooled_Bac_headers.sam
  rm ${MD_WORKDIR_P2_OUTDIR}/Pooled_bac.sam
done

echo "MD Phase IIb: Combine sam with the headers"
cat ${MD_WORKDIR_P2_OUTDIR}/Pooled_Bac_headers.sam \
    ${MD_WORKDIR_P2_OUTDIR}/Pooled_Bac_noHeaders.sam \
    > ${MD_WORKDIR_P2_OUTDIR}/Pooled_Bac.sam

rm ${MD_WORKDIR_P2_OUTDIR}/Pooled_Bac_headers.sam
rm ${MD_WORKDIR_P2_OUTDIR}/Pooled_Bac_noHeaders.sam

echo "MD Phase IId: Look for best hits between reads in updated sam file"

python ${pathoscope_python} \
-t sam -e csmd -f ${MD_WORKDIR_P2_OUTDIR}/Pooled_Bac.sam \
-outdir ${MD_WORKDIR_P2_OUTDIR}

rm ${MD_WORKDIR_P2_OUTDIR}/Pooled_Bac.sam
