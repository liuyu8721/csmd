# All rights reserved by Liuyu
# Author: Liuyu
#!/bin/sh

do_download_AHG () {
  db_path=$1
  if [ -d "${db_path}" ]
  then
    ahg_dir=${db_path}/AHG
    ahg_download_dir=${ahg_dir}/SEQ

    if [ -d "${ahg_download_dir}" ]; then
     echo -e "\033[44;37;5m ERROR:\033[0m path ${ahg_download_dir} has already existed!"
     exit
    fi

    mkdir -p ${ahg_download_dir}
    cd ${ahg_download_dir}
    echo "downloading HuRef"
    wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/002/125/GCA_000002125.2_HuRef/GCA_000002125.2_HuRef_genomic.fna.gz -P ${ahg_download_dir}
    echo "downloading YH"
    wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/004/845/GCA_000004845.2_YH_2.0/GCA_000004845.2_YH_2.0_genomic.fna.gz -P ${ahg_download_dir}
    echo "downloading BGIYH"
    wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/465/GCA_000005465.1_BGIAF/GCA_000005465.1_BGIAF_genomic.fna.gz -P ${ahg_download_dir}

    echo "concatenating all fna files into ahg.fasta"
    zcat *.fna.gz > ${ahg_dir}/ahg.fasta
  else
    echo "please provide db path and make sure its path exist"
  fi
}

do_download_EHG () {
  db_path=$1
  if [ -d "${db_path}" ]
  then
    ehg_dir=${db_path}/EHG
    ehg_download_dir=${ehg_dir}/SEQ

    if [ -d "${ehg_download_dir}" ]; then
     echo -e "\033[44;37;5m ERROR:\033[0m path ${ehg_download_dir} has already existed!"
     exit
    fi

    mkdir -p ${ehg_download_dir}
    cd ${ehg_download_dir}
    echo "downloading 'Ensembl Homo sapiens cDNA database'"
    wget -c ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz \
      -P ${ehg_download_dir}
    echo "downloading 'NCBI Homo sapiens RNA database'"
    wget -c ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/RNA/rna.fa.gz \
      -P ${ehg_download_dir}
    echo "downloading 'NCBI BLAST human genome database'"
    ### Downloading the Blast index files
    wget -c ftp://ftp.ncbi.nlm.nih.gov/blast/db/human_genomic.*.tar.gz \
      -P ${ehg_download_dir}
    ### Extract fasta from the index files
    for item in `ls ${ehg_download_dir}/human_genomic.*.tar.gz`;
    do
        tar xzvf ${item}
    done
    blastdbcmd -entry all -db ${ehg_download_dir}/human_genomic -out ${ehg_download_dir}/NCBIBlastHumanSequence.fasta
    rm ${ehg_download_dir}/human_genomic.??.n* 
    rm ${ehg_download_dir}/human_genomic.nal
    ### Split the fasta file in order that each one has the size less than 6G.
    grep '^>' ${ehg_download_dir}/NCBIBlastHumanSequence.fasta | sed 's/.//' > ${ehg_download_dir}/header.list
    sed -n '1,128p' ${ehg_download_dir}/header.list > ${ehg_download_dir}/ehg.split.list00
    sed -n '129,558p' ${ehg_download_dir}/header.list > ${ehg_download_dir}/ehg.split.list01
    sed -n '559,846p' ${ehg_download_dir}/header.list > ${ehg_download_dir}/ehg.split.list02
    sed -n '847,1294p' ${ehg_download_dir}/header.list > ${ehg_download_dir}/ehg.split.list03
    sed -n '1295,1683p' ${ehg_download_dir}/header.list > ${ehg_download_dir}/ehg.split.list04
    sed -n '1684,1835p' ${ehg_download_dir}/header.list > ${ehg_download_dir}/ehg.split.list05
    sed -n '1836,1998p' ${ehg_download_dir}/header.list > ${ehg_download_dir}/ehg.split.list06
    sed -n '1999,2171p' ${ehg_download_dir}/header.list > ${ehg_download_dir}/ehg.split.list07
    sed -n '2172,2421p' ${ehg_download_dir}/header.list > ${ehg_download_dir}/ehg.split.list08
    sed -n '2422,3472p' ${ehg_download_dir}/header.list > ${ehg_download_dir}/ehg.split.list09
    sed -n '3473,3505p' ${ehg_download_dir}/header.list > ${ehg_download_dir}/ehg.split.list10
    for item in `ls ${ehg_download_dir}/ehg.split.list*`
    do
        no=$(echo ${item} | awk -F "." '{print $NF}' | cut -c5-6)
        seqtk subseq -l 80 ${ehg_download_dir}/NCBIBlastHumanSequence.fasta ${item} \
          > ${ehg_dir}/ehg.${no}.fasta
    done
    zcat ${ehg_download_dir}/Homo_sapiens.GRCh38.cdna.all.fa.gz ${ehg_download_dir}/rna.fa.gz \
        > ${ehg_dir}/ehg.11.fasta
    rm ${ehg_download_dir}/human_genomic.*
    gzip ${ehg_download_dir}/NCBIBlastHumanSequence.fasta
  else
    echo "please provide db path and make sure its path exist"
  fi
}

do_download_nt () {
  db_path=$1
  if [ -d "${db_path}" ]
  then
    nt_dir=${db_path}/nt
    nt_download_dir=${db_path}/nt/SEQ

    if [ -d "${nt_download_dir}" ]; then
     echo -e "\033[44;37;5m ERROR:\033[0m path ${nt_download_dir} has already existed!"
     exit
    fi

    mkdir -p ${nt_download_dir}
    cd ${nt_download_dir}
    wget -c ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz \
      -P ${nt_download_dir}
    gzip -d ${nt_download_dir}/nt.gz
  else
    echo "please provide db path and make sure its path exist"
  fi
}

do_download_RefSeqBac () {
  db_path=$1
  if [ -d "${db_path}" ]
  then
    bacteria_dir=${db_path}/RefSeq/bacteria

    if [ -d "${bacteria_dir}" ]; then
     echo -e "\033[44;37;5m ERROR:\033[0m path ${bacteria_dir} has already existed!"
     exit
    fi

    mkdir -p ${bacteria_dir}
    sh ${CSMD_HOME}/MD_P0.sh ${bacteria_dir}
  else
    echo "please provide db path and make sure its path exist"
  fi
}

do_download_taxonomy () {
  db_path=$1
  if [ -d "${db_path}" ]
  then
    # Edited by Liu
    taxonomy_dir=${db_path}/taxonomy/taxdump

    if [ -d "${taxonomy_dir}" ]; then
     echo -e "\033[44;37;5m ERROR:\033[0m path ${taxonomy_dir} has already existed!"
     exit
    fi

    mkdir -p ${taxonomy_dir}
    cd ${taxonomy_dir}
    wget -c ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz \
      -P ${taxonomy_dir}
    tar -xvf taxdump.tar.gz -C ${taxonomy_dir}
  else
    echo "please provide db path and make sure its path exist"
  fi
}

do_download_hg38 () {
  db_path=$1
  if [ -d "${db_path}" ]
  then
    hg_38_dir=${db_path}/hg38
    hg_38_download_dir=${hg_38_dir}/SEQ

    if [ -d "${hg_38_download_dir}" ]; then
     echo -e "\033[44;37;5m ERROR:\033[0m path ${hg_38_download_dir} has already existed!"
     exit
    fi

    mkdir -p ${hg_38_download_dir}
    cd ${hg_38_download_dir}
    wget -c ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/*  -P ${hg_38_download_dir}
    for item in `ls ${hg_38_download_dir}/*.fa.gz`;
    do
      gunzip ${item}
    done
    echo "concatenating all fa files into hg38_ucsc.fasta"
    cat *.fa > ${hg_38_dir}/hg38_ucsc.fasta
  else
    echo "please provide db path and make sure its path exist"
  fi
}
