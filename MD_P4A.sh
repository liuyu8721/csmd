# All rights reserved by Liuyu
# Author: Liuyu
#!/bin/sh
source ${CSMD_HOME}/configs.sh

MD_WORKDIR_P4A_SEQDIR=$1
MD_WORKDIR_P4A_SEQLST=$2
MD_WORKDIR_P4A_OUTDIR=$3

echo "#--------------------------------------------------------------#"
#echo "MD Phase IVa: Genome collection"

for aa in $(cat ${MD_WORKDIR_P4A_SEQLST})
do
    taxid=$(awk -F '\t' '{if($1=="'${aa}'")print $6}' ${MD_WORKDIR_P4A_SEQDIR}/assembly_summary.txt)
    organism=$(awk -F '\t' '{if($1=="'${aa}'")print $8}' ${MD_WORKDIR_P4A_SEQDIR}/assembly_summary.txt | sed 's/ /_/g' | sed 's/\//!/g')
    header=">${aa}|${taxid}|${organism}"
    echo ${header}

    fasta_genome=$(ls ${MD_WORKDIR_P4A_SEQDIR}/SEQ/${aa}*.gz)
    if [ "$fasta_genome" != "" ]; then
	gzip -d ${fasta_genome} -c > ${MD_WORKDIR_P4A_OUTDIR}/DB/Genome/${aa}.fna
	sed -i 's/^>.*$/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/' ${MD_WORKDIR_P4A_OUTDIR}/DB/Genome/${aa}.fna
	sed -i "1s/^.*$/${header}/" ${MD_WORKDIR_P4A_OUTDIR}/DB/Genome/${aa}.fna
    else
	echo -e "\033[44;37;5m ERROR:\033[0m the fasta genome of ${aa} not exist"
	exit
    fi
done

cat ${MD_WORKDIR_P4A_OUTDIR}/DB/Genome/*.fna > ${MD_WORKDIR_P4A_OUTDIR}/DB/CorrectSpecies
cd ${MD_WORKDIR_P4A_OUTDIR}/DB
bowtie2-build CorrectSpecies CorrectSpecies;

echo "#--------------------------------------------------------------#"
