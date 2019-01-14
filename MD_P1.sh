# All rights reserved by Liuyu
# Author: Liuyu
#!/bin/sh
source ${CSMD_HOME}/configs.sh

MD_WORKDIR_P1=$1

cd `dirname $0`

echo "#--------------------------------------------------------------#"
echo "MD Phase I: Generate an initial redundant reference library"
if [ -f "${MD_WORKDIR_P1}/assembly_summary.txt" ];then
  Rscript ${CSMD_HOME}/RepSpecies.R ${MD_WORKDIR_P1}
  Rscript ${CSMD_HOME}/RepGenomes.R ${MD_WORKDIR_P1}
else 
  echo -e "\033[44;37;5m ERROR: \033[0m ${MD_WORKDIR_P1}/assembly_summary.txt not exit"
  exit
fi

# It will takes some time to finish the job
if [ ! -d "${MD_WORKDIR_P1}/SEQ/" ];then
  echo -e "\033[44;37;5m ERROR:\033[0m ${MD_WORKDIR_P1}/SEQ/ not exist"
  echo "Please make sure all RefSeq bacteria genomes can be found in ${MD_WORKDIR_P1}/SEQ"
  exit
fi

if [ ! -d "${MD_WORKDIR_P1}/RepSpecies/Genome/" ];then
  mkdir -p ${MD_WORKDIR_P1}/RepSpecies/Genome
else
  if [ "`ls -A ${MD_WORKDIR_P1}/RepSpecies/Genome/`" != "" ];then
	  echo -e "\033[44;37;5m ERROR:\033[0m ${MD_WORKDIR_P1}/RepSpecies/Genome/ not empty"
	  exit
  fi
fi

cd ${MD_WORKDIR_P1}
grep "^GCF" SpeciesRepresentative.csv | \
    split -l 6 -d -a 2 - SpeciesRepresentative.split.
for RepSpecies in $(ls SpeciesRepresentative.split.??)
do
  seqno=$(echo ${RepSpecies} | cut -d"." -f3)
  mkdir -p ${MD_WORKDIR_P1}/RepSpecies/Genome/${seqno}
  for aa in $(cat ${RepSpecies} | awk -F"," '{print $1}')
  do
	  taxid=$(cat ${RepSpecies} | grep "^${aa}" | awk -F"," '{print $6}')
	  organism=$(cat ${RepSpecies} | grep "^${aa}" | awk -F"," '{print $8}' | sed 's/ /_/g' | sed 's/\//!/g')
	  header=">${aa}|${taxid}|${organism}"
	  echo ${header}
	fasta_genome=$(ls ${MD_WORKDIR_P1}/SEQ/10000ok/${aa}*.gz)
	if [ "$fasta_genome" != "" ]; then
		gzip -d ${fasta_genome} -c > ${MD_WORKDIR_P1}/RepSpecies/Genome/${seqno}/${aa}.fna
	else
		echo -e "\033[44;37;5m ERROR:\033[0m the fasta genome of ${aa} not exit"
		exit
	fi
	sed -i 's/^>.*$/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/'  ${MD_WORKDIR_P1}/RepSpecies/Genome/${seqno}/${aa}.fna
	sed -i "1s/^.*$/${header}/" ${MD_WORKDIR_P1}/RepSpecies/Genome/${seqno}/${aa}.fna
  done
  cat ${MD_WORKDIR_P1}/RepSpecies/Genome/${seqno}/*.fna > ${MD_WORKDIR_P1}/RepSpecies/RefSeq_RepSpecies_${seqno}
  cd ${MD_WORKDIR_P1}/RepSpecies
  bowtie2-build RefSeq_RepSpecies_${seqno} RefSeq_RepSpecies_${seqno}
  rm RefSeq_RepSpecies_${seqno}
  cd ${MD_WORKDIR_P1}
  rm ${RepSpecies}
done

if [ ! -d "${MD_WORKDIR_P1}/RepGenomes/Genome/" ];then
  mkdir -p ${MD_WORKDIR_P1}/RepGenomes/Genome
else
  if [ "`ls -A ${MD_WORKDIR_P1}/RepGenomes/Genome/`" != "" ];then
    echo -e "\033[44;37;5m ERROR:\033[0m ${MD_WORKDIR_P1}/RepGenomes/Genome/ not empty"
    exit
  fi
fi

for aa in $(cat GenomeRepresentative.csv | grep "^GCF" | awk -F"," '{print $1}')
do
	fasta_genome=$(ls ${MD_WORKDIR_P1}/SEQ/10000ok/${aa}*.gz)
    if [ "$fasta_genome" != "" ]; then
      gzip -d ${fasta_genome} -c > ${MD_WORKDIR_P1}/RepGenomes/Genome/${aa}.fna
    else
      echo -e "\033[44;37;5m ERROR:\033[0m the fasta genome of ${aa} not exit"
      exit
    fi
#	gzip -d ${MD_WORKDIR_P1}/SEQ/10000ok/${aa}*.gz -c > ${MD_WORKDIR_P1}/RepGenomes/Genome/${aa}.fna
done

cat ${MD_WORKDIR_P1}/RepGenomes/Genome/*.fna > ${MD_WORKDIR_P1}/RepGenomes/RefSeq_RepGenomes
cd ${MD_WORKDIR_P1}/RepGenomes
makeblastdb -in RefSeq_RepGenomes -dbtype nucl
rm RefSeq_RepGenomes
