# All rights reserved by Liuyu
# Author: Liuyu
#!/bin/sh
source ${CSMD_HOME}/configs.sh

MD_WORKDIR_P3_CUTOFF=$1
MD_WORKDIR_P3_THREAD=$2
MD_WORKDIR_P3_SAM=$3

MD_WORKDIR_P3_REPORT=$4

MD_WORKDIR_P3_NTNAME=$5

MD_WORKDIR_P3_RSNAME=$6

MD_WORKDIR_P3_TAXDIR=$7
MD_WORKDIR_P3_OUTDIR=$8
echo "#--------------------------------------------------------------#"
echo "MD Phase III: Species correction"
if [ ! -f "${MD_WORKDIR_P3_SAM}" ];then
  echo -e "\033[44;37;5m ERROR: \033[0m the csmd discovery sam file ${MD_WORKDIR_P3_SAM} not exist"
  exit
fi

if [ ! -f "${MD_WORKDIR_P3_REPORT}" ];then
  echo -e "\033[44;37;5m ERROR: \033[0m the csmd discovery report file ${MD_WORKDIR_P3_REPORT} not exit"
  exit
fi

if [ ! -d "${MD_WORKDIR_P3_OUTDIR}" ];then
  echo -e "\033[44;37;5m ERROR:\033[0m OUTPATH can not be found"
  exit
elif [ "`ls -A ${MD_WORKDIR_P3_OUTDIR}`" != "" ];then
  echo -e "\033[44;37;5m ERROR:\033[0m ${MD_WORKDIR_P3_OUTDIR} not empty"
  exit
fi

if [[ ! -f "${MD_WORKDIR_P3_TAXDIR}/taxdb.btd" ]] || [[ ! -f "${MD_WORKDIR_P3_TAXDIR}/taxdb.bti" ]];then
  echo -e "\033[44;37;5m ERROR: \033[0m taxdb.btd or taxdb.bti not exist in the ${MD_WORKDIR_P3_TAXDIR}"
  exit
fi

if [[ ! -f "${MD_WORKDIR_P3_TAXDIR}/taxtree.txt" ]] || [[ ! -f "${MD_WORKDIR_P3_TAXDIR}/GenBank2RefSeq.txt" ]];then
  echo -e "\033[44;37;5m ERROR: \033[0m taxtree.txt or GenBank2RefSeq.txt exist in the ${MD_WORKDIR_P3_TAXDIR}"
  exit
fi

cd ${MD_WORKDIR_P3_OUTDIR}

echo "MD Phase IIIa: Ready files for species correction"
samtools view -S -H ${MD_WORKDIR_P3_SAM} > csmd_discovery_headers.sam
grep "^GCF_" ${MD_WORKDIR_P3_REPORT} | \
awk -F '\t' '{if($4>'"${MD_WORKDIR_P3_CUTOFF}"')print $1;}' | cut -d'|' -f1 \
> screening_species.list

echo "MD Phase IIIb: Blast analysis"
cp ${MD_WORKDIR_P3_TAXDIR}/taxdb.btd .
cp ${MD_WORKDIR_P3_TAXDIR}/taxdb.bti .
mkdir ${MD_WORKDIR_P3_OUTDIR}/FASTQ
mkdir ${MD_WORKDIR_P3_OUTDIR}/FASTA
mkdir ${MD_WORKDIR_P3_OUTDIR}/BLAST
mkdir ${MD_WORKDIR_P3_OUTDIR}/BLAST/nt
mkdir ${MD_WORKDIR_P3_OUTDIR}/BLAST/RefSeq
cat screening_species.list | while read line; do

	id=${line}
	echo ${id}

	### Extract mapped reads
	samtools view -S ${MD_WORKDIR_P3_SAM} | grep "${id}" > ${id}_mappedReads.sam
	cat csmd_discovery_headers.sam ${id}_mappedReads.sam > ${id}.sam
	samtools view -S -b ${id}.sam -o ${id}.bam

	rm ${id}_mappedReads.sam
	rm ${id}.sam

	### Bam to fastq
	bamToFastq -i ${id}.bam -fq ${id}.fastq
	rm ${id}.bam

	### Sampling
	fastq-sample -n 1000 -o ${id}.sample ${id}.fastq

	### Fastq to Fasta
	module load seqtk/master
	seqtk seq -A ${id}.sample.fastq > ${id}.sample.fasta
	rm ${id}.sample.fastq

	### Blast based on nt database
	blastn \
	-db ${MD_WORKDIR_P3_NTNAME} \
	-query ${id}.sample.fasta \
	-num_threads ${MD_WORKDIR_P3_THREAD} -task megablast \
	-outfmt '6 std staxids' -out ${id}.sample.nt.blast \
	-qcov_hsp_perc 50 -max_target_seqs 100 -perc_identity 90 -evalue 1e-5

	### Blast based on RefSeq database
	blastn \
	-db ${MD_WORKDIR_P3_RSNAME} \
	-query ${id}.sample.fasta \
	-num_threads ${MD_WORKDIR_P3_THREAD} -task megablast \
	-outfmt 6 -out ${id}.sample.RefSeq.blast \
	-qcov_hsp_perc 50 -max_target_seqs 100 -perc_identity 90 -evalue 1e-5

	mv ${id}.fastq ${MD_WORKDIR_P3_OUTDIR}/FASTQ
	mv ${id}.sample.fasta ${MD_WORKDIR_P3_OUTDIR}/FASTA
	mv ${id}.sample.nt.blast ${MD_WORKDIR_P3_OUTDIR}/BLAST/nt
	mv ${id}.sample.RefSeq.blast ${MD_WORKDIR_P3_OUTDIR}/BLAST/RefSeq
done
rm taxdb.btd taxdb.bti csmd_discovery_headers.sam

echo "MD Phase IIIc: Summarize the blast results"
Rscript ${CSMD_HOME}/BlastAnalysisForNT.R ${MD_WORKDIR_P3_TAXDIR} ${MD_WORKDIR_P3_OUTDIR} ${MD_WORKDIR_P3_CUTOFF}
Rscript ${CSMD_HOME}/BlastAnalysisForRefSeq.R ${MD_WORKDIR_P3_TAXDIR} ${MD_WORKDIR_P3_OUTDIR} ${MD_WORKDIR_P3_CUTOFF}
Rscript ${CSMD_HOME}/BlastAnalysisSummary.R ${MD_WORKDIR_P3_TAXDIR} ${MD_WORKDIR_P3_OUTDIR} ${MD_WORKDIR_P3_CUTOFF}
