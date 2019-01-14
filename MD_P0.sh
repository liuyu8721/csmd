# All rights reserved by Liuyu
# Author: Liuyu
#
# This shell script is to download all genomic fna files to local directory.
# e.g. sh ./MD_P0.sh /data/dna/RefSeq/bacteria
#!/bin/bash
source ${CSMD_HOME}/configs.sh

#REF_SEQ_DB_PATH=/dbpath/RefSeq/bacteria
REF_SEQ_DB_PATH=$1
REF_SEQ_DB_DOWNLOAD_PATH=${REF_SEQ_DB_PATH}/SEQ

remind () {
  echo "Please be patient since this takes quite long time to finish."
  echo "MD Phase Ia: Collect and modify the FTP URLs to point to the _genomic.fna.gz files"
}

prepare_dir () {
	if [ -d "${REF_SEQ_DB_DOWNLOAD_PATH}" ]; then
	  echo "${REF_SEQ_DB_DOWNLOAD_PATH} already exists, please check to delete this dir and run again"
    exit;
  fi
}

download_summary_file () {
	echo "Download assembly_summary.txt to prepare for downloading data"
  mkdir -p ${REF_SEQ_DB_DOWNLOAD_PATH}
  cd ${REF_SEQ_DB_DOWNLOAD_PATH}
  wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt -P ${REF_SEQ_DB_PATH}
}

collect_file_list () {
  awk '{FS="\t"} !/^#/ {print $20} ' ${REF_SEQ_DB_PATH}/assembly_summary.txt | \
  sed -r 's|(ftp://ftp.ncbi.nlm.nih.gov/genomes/all/.+/)(GCF_.+)|\1\2/\2_genomic.fna.gz|' > genomic_file.txt
  echo "Done with collecting data file list."	
}

remove_old_data_files () {
	echo "remove all old data files."
	rm -rf *.fna.gz
	rm -rf *.fna
}

download_data_files () {
	echo "start to download all data files."
	for item in $(cat genomic_file.txt); do wget -c ${item}; done
}

unzip_files () {
	gunzip *.fna.gz
}

do_proceed () {
  remind
  prepare_dir
  download_summary_file
  cd ${REF_SEQ_DB_DOWNLOAD_PATH}
  collect_file_list
  remove_old_data_files
  download_data_files
  unzip_files
}

execute_md_p0 () {
  echo "MD Phase 0: will download all _genomic.fna.gz files from remote FTP."
  echo "Please confirm that ${REF_SEQ_DB_DOWNLOAD_PATH} data will be erased!!!"
  read -p "Continue (y/n)?" choice
  case "$choice" in 
    y|Y ) do_proceed;;
    n|N ) echo "no";;
    * ) echo "invalid";;
  esac
}

execute_md_p0





