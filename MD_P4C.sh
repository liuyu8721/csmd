# All rights reserved by Liuyu
# Author: Liuyu
#!/bin/sh
source ${CSMD_HOME}/configs.sh

MD_WORKDIR_P4C_INPUT=$1
MD_WORKDIR_P4C_OUTDIR=$2
MD_WORKDIR_P4C_CODEDIR=$3

echo "#--------------------------------------------------------------#"
#echo "MD Phase IVc: Coverage structure analysis"
cd ${MD_WORKDIR_P4C_OUTDIR}/COV

samtools view -S ${MD_WORKDIR_P4C_INPUT} > mappedReads.sam
samtools view -S -H ${MD_WORKDIR_P4C_INPUT} > header.sam

grep "^@SQ" header.sam | awk '{print $2}' | cut -d"|" -f1 | cut -d":" -f2 > rsid
grep "^@SQ" header.sam | awk '{print $2}' | cut -d"|" -f3 > taxname
grep "^@SQ" header.sam | awk '{print $3}' | cut -d":" -f2 > len
paste rsid taxname len > speciesList.txt
rm rsid taxname len

cat speciesList.txt | while read line; do 
	species=$(echo $line | awk '{print $2}')
	speShort=$(echo $line | awk '{print $1}')
	genomeLength=$(echo $line | awk '{print $3}')
	echo ${speShort} ${species} ${genomeLength}  

	grep "${speShort}" mappedReads.sam > temp.sam
	cat header.sam temp.sam > temp1.sam
	samtools view -S -b temp1.sam -o temp1.bam
	samtools sort temp1.bam -o temp2.bam
	samtools index temp2.bam
	samtools view temp2.bam -o ${speShort}
	Rscript ${MD_WORKDIR_P4C_CODEDIR}/covBins.R ${speShort} ${genomeLength} 5000 
	rm temp.sam temp1.sam temp1.bam temp2.bam ${speShort} temp2.bam.bai

     	grep "${speShort}" mappedReads.sam | grep -v "XS" > temp.sam
     	cat header.sam temp.sam > temp1.sam
     	samtools view -S -b temp1.sam -o temp1.bam
     	samtools sort temp1.bam -o temp2.bam
     	samtools index temp2.bam
     	samtools view temp2.bam -o ${speShort}.unique
     	Rscript ${MD_WORKDIR_P4C_CODEDIR}/covBins.R ${speShort}.unique ${genomeLength} 5000
    	rm temp.sam temp1.sam temp1.bam temp2.bam ${speShort}.unique temp2.bam.bai
done
rm mappedReads.sam header.sam speciesList.txt

echo "#--------------------------------------------------------------#"

