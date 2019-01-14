# All rights reserved by Liuyu
# Author: Liuyu
#!/bin/sh
source ${CSMD_HOME}/configs.sh

MD_WORKDIR_P4_THREAD=$1
MD_WORKDIR_P4_SEQDIR=$2
MD_WORKDIR_P4_SEQLST=$3
MD_WORKDIR_P4_INPUT=$4
MD_WORKDIR_P4_OUTDIR=$5

echo "${MD_WORKDIR_P4_OUTDIR}"
echo "#--------------------------------------------------------------#"
echo "MD Phase IV: Species refinement"
if [ ! -d "${MD_WORKDIR_P4_SEQDIR}/SEQ/" ];then
    echo -e "\033[44;37;5m ERROR:\033[0m ${MD_WORKDIR_P4_SEQDIR}/SEQ/ not exist"
    echo "Please make sure all RefSeq bacteria genomes can be found in ${MD_WORKDIR_P4_SEQDIR}/SEQ/"
    exit
elif [ ! -f "${MD_WORKDIR_P4_SEQDIR}/assembly_summary.txt" ];then
    echo -e "\033[44;37;5m ERROR: \033[0m ${MD_WORKDIR_P4_SEQDIR}/assembly_summary.txt not exit"
    exit
else
    echo -e "\033[44;37;5m TIPS:\033[0m The SEQDIR is ready."
    echo "NOTE: All RefSeq bacteria genomes are expected in ${MD_WORKDIR_P4_SEQDIR}/SEQ/"
    echo "NOTE: RefSeq bacteria summary information is expected in ${MD_WORKDIR_P4_SEQDIR}/assembly_summary.txt"
fi

if [ ! -f "${MD_WORKDIR_P4_SEQLST}" ];then
    echo -e "\033[44;37;5m ERROR:\033[0m the sequence list ${MD_WORKDIR_P4_SEQLST} not exist"
    echo "Please make sure the sequence list is ready for database update, each line with a RefSeq accssion no. "
else
    echo -e "\033[44;37;5m TIPS:\033[0m The SEQLST is ready."
    echo "NOTE: Sequences for database update are expected in ${MD_WORKDIR_P4_SEQLST}, each line with a RefSeq accssion no."
fi

if [ ! -f "${MD_WORKDIR_P4_INPUT}" ];then
    echo -e "\033[44;37;5m ERROR: \033[0m ${MD_WORKDIR_P4_INPUT} not exist"
    exit
fi

echo "MD Phase IVa: Genome collection"
if [ ! -d "${MD_WORKDIR_P4_OUTDIR}/DB/Genome/" ];then
    mkdir -p ${MD_WORKDIR_P4_OUTDIR}/DB/Genome
else
    if [ "`ls -A ${MD_WORKDIR_P4_OUTDIR}/DB/Genome/`" != "" ];then
        echo -e "\033[44;37;5m ERROR:\033[0m ${MD_WORKDIR_P4_OUTDIR}/DB/Genome/ not empty"
        exit
    fi
fi

${CSMD_HOME}/MD_P4A.sh ${MD_WORKDIR_P4_SEQDIR} ${MD_WORKDIR_P4_SEQLST} ${MD_WORKDIR_P4_OUTDIR}

echo "MD Phase IVb: Re-alignment to the udpated database"
bowtie2 \
-p ${MD_WORKDIR_P4_THREAD} -k 10 --very-sensitive \
-x ${MD_WORKDIR_P4_OUTDIR}/DB/CorrectSpecies \
-f ${MD_WORKDIR_P4_INPUT} \
-S ${MD_WORKDIR_P4_OUTDIR}/Pooled_bac_updateDB.sam

echo "MD Phase IVc: Coverage structure analysis"
if [ ! -d "${MD_WORKDIR_P4_OUTDIR}/COV/" ];then
    mkdir -p ${MD_WORKDIR_P4_OUTDIR}/COV
else
    if [ "`ls -A ${MD_WORKDIR_P4_OUTDIR}/COV/`" != "" ];then
        echo -e "\033[44;37;5m ERROR:\033[0m ${MD_WORKDIR_P4_OUTDIR}/COV/ not empty"
       exit
    fi
fi

if [ ! -f "${MD_WORKDIR_P4_OUTDIR}/Pooled_bac_updateDB.sam" ];then
    echo -e "\033[44;37;5m ERROR: \033[0m ${MD_WORKDIR_P4_OUTDIR}/Pooled_bac_updateDB.sam not ready"
    exit
fi

${CSMD_HOME}/MD_P4C.sh ${MD_WORKDIR_P4_OUTDIR}/Pooled_bac_updateDB.sam ${MD_WORKDIR_P4_OUTDIR} ${CSMD_HOME}

Rscript ${CSMD_HOME}/CovAnalysis.R ${CSMD_HOME} ${MD_WORKDIR_P4_OUTDIR}/COV
