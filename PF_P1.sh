# All rights reserved by Liuyu
# Author: Liuyu
#!/bin/sh
source ${CSMD_HOME}/configs.sh

PF_WORKDIR_P1_SEQDIR=$1
PF_WORKDIR_P1_SEQLST=$2
PF_WORKDIR_P1_OUTDIR=$3

echo "#--------------------------------------------------------------#"
echo "Per-sample microbiome profiling"

if [ ! -d "${PF_WORKDIR_P1_SEQDIR}/SEQ/" ];then
    echo -e "\033[44;37;5m ERROR:\033[0m ${PF_WORKDIR_P1_SEQDIR}/SEQ/ not exist"
    echo "Please make sure all RefSeq bacteria genomes can be found in ${PF_WORKDIR_P1_SEQDIR}/SEQ/"
    exit
elif [ ! -f "${PF_WORKDIR_P1_SEQDIR}/assembly_summary.txt" ];then
    echo -e "\033[44;37;5m ERROR: \033[0m ${PF_WORKDIR_P1_SEQDIR}/assembly_summary.txt not exit"
    exit
else
    echo -e "\033[44;37;5m TIPS:\033[0m The SEQDIR is ready."
    echo "NOTE: All RefSeq bacteria genomes are expected in ${PF_WORKDIR_P1_SEQDIR}/SEQ/"
    echo "NOTE: RefSeq bacteria summary information is expected in ${PF_WORKDIR_P1_SEQDIR}/assembly_summary.txt"
fi

if [ ! -f "${PF_WORKDIR_P1_SEQLST}" ];then
    echo -e "\033[44;37;5m ERROR:\033[0m the sequence list ${PF_WORKDIR_P1_SEQLST} not exist"
    echo "Please make sure the sequence list is ready for database update, each line with a RefSeq accssion no. "
else
    echo -e "\033[44;37;5m TIPS:\033[0m The SEQLST is ready."
    echo "NOTE: Sequences for database update are expected in ${PF_WORKDIR_P1_SEQLST}, each line with a RefSeq accssion no."
fi

echo "PF Phase I: Genome colection and index"
if [ ! -d "${PF_WORKDIR_P1_OUTDIR}/DB/Genome/" ];then
   mkdir -p ${PF_WORKDIR_P1_OUTDIR}/DB/Genome
else
   if [ "`ls -A ${PF_WORKDIR_P1_OUTDIR}/DB/Genome/`" != "" ];then
       echo -e "\033[44;37;5m ERROR:\033[0m ${PF_WORKDIR_P1_OUTDIR}/DB/Genome/ not empty"
       exit
   fi
fi

for aa in $(cat ${PF_WORKDIR_P1_SEQLST})
do
    taxid=$(awk -F '\t' '{if($1=="'${aa}'")print $6}' ${PF_WORKDIR_P1_SEQDIR}/assembly_summary.txt)
    organism=$(awk -F '\t' '{if($1=="'${aa}'")print $8}' ${PF_WORKDIR_P1_SEQDIR}/assembly_summary.txt | sed 's/ /_/g' | sed 's/\//!/g')
    header=">${aa}|${taxid}|${organism}"
    echo ${header}

    fasta_genome=$(ls ${PF_WORKDIR_P1_SEQDIR}/SEQ/${aa}*.gz)
    if [ "$fasta_genome" != "" ]; then
	gzip -d ${fasta_genome} -c > ${PF_WORKDIR_P1_OUTDIR}/DB/Genome/${aa}.fna
       	sed -i 's/^>.*$/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/' ${PF_WORKDIR_P1_OUTDIR}/DB/Genome/${aa}.fna
	sed -i "1s/^.*$/${header}/" ${PF_WORKDIR_P1_OUTDIR}/DB/Genome/${aa}.fna
    else
	echo -e "\033[44;37;5m ERROR:\033[0m the fasta genome of ${aa} not exist"
	exit
    fi
done

cat ${PF_WORKDIR_P1_OUTDIR}/DB/Genome/*.fna > ${PF_WORKDIR_P1_OUTDIR}/DB/csmdSpecies
cd ${PF_WORKDIR_P1_OUTDIR}/DB
bowtie2-build csmdSpecies csmdSpecies;
