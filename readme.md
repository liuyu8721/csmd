Computation Subtraction-based Microbiome Discovery (CSMD)
====
CSMD is a computational pipeline for high-resolution profiling of low abundance microbiome in clinical samples using whole genome shotgun sequencing.
## 1. Pre-requisite software
Table 1. List of pre-requisite software and the available information

|No.            |Software       |Version        | Availability |
| ------------- |-------------|-------------| -----       |
|1|BWA|0.7.15|http://bio-bwa.sourceforge.net/|
|2|SAMtools|1.4|http://www.htslib.org/|
|3|bedtools|2.26.0|https://bedtools.readthedocs.io/en/latest/|
|4|seqtk|1.2|https://github.com/lh3/seqtk|
|5|RepeatMasker|4.0.7|http://repeatmasker.org/|
|6|PathoScope 2.0|0.02|https://sourceforge.net/projects/pathoscope/|
|7|R|3.3.3|https://www.r-project.org/|
|8|countreg|0.2-0|https://r-forge.r-project.org/projects/countreg/|
|9|Bowtie2|2.2.1|https://sourceforge.net/projects/bowtie-bio/files/bowtie2/|
|10|fastq-tools|0.8|https://github.com/dcjones/fastq-tools|
|11|blast|2.6.0+|https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download|
## 2. Database availability
Table 2. List of pre-download databases and the available information

|No.|Databases|Availability|
|----|----|----|
|1|Human reference genome (hg38)|https://genome.ucsc.edu/cgi-bin/hgGateway?db=hg38|
|2|Three assembled human genomes available on NCBI: HuRef, YH and BGIAF.|ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/002/125/GCA_000002125.2_HuRef<br/> ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/004/845/GCA_000004845.2_YH_2.0<br/> ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/465/GCA_000005465.1_BGIAF</br>|
|3|Repbase|https://www.girinst.org/|
|4|Ensembl Homo sapiens cDNA database|ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/|
|5|NCBI Homo sapiens RNA database|ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/RNA/|
|6|NCBI BLAST human genome database|ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/|
|7|NCBI non-redundant nucleotide sequences (nt)|ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/|
|8|NCBI RefSeq Bacteria database|ftp://ftp.ncbi.hlm.nih.gov/genomes/refseq/bacteria|
|9|Taxonomy files|https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz|
## 3. Computational pipeline
### Step 3.1: A four-phase human-derived sequence subtraction.
* Phase I: Subtract reads from standard human reference genome (hg38)
* Phase II: Subtract reads from three additional assembled human genomes (AHG), including HuRef, YH and BGIAF
* Phase III: Subtract low complexity reads (LCR)
* Phase IV: Subtract reads from three extra human sequence databases (EHG), including Ensembl Homo sapiens cDNA database, NCBI Homo sapiens RNA database, and NCBI BLAST human genome database
### Step 3.2: A four-phase microbiome discovery procedure.
* Phase I: Generate an initial redundant reference library that may include all the species genomes of interest.
* Phase II: Identify a list of alternatives of candidate genomes through fast similarity search against the initial library.
* Phase III: Screen out genomes with significantly insufficient coverage and do the species correction for possibly misidentified genomes using BLAST analysis.
* Phase IV: Perform genome refinement through analyzing their coverage structure.
### Step 3.3: Per-sample microbiome profiling.
## 4.Installation
### 4.1 Download
Download the code from https://github.com/liuyu8721/CSMD. 
You could issue the following command to extract the files: "unzip csmd-master.zip".

### 4.2 Configuration
1. Make sure the pre-requisite software listed in Table 1 has been installed and available in the runtime environment. They should be first configured in 'config.sh', for example:
```shell
export PATH=/public/software/bwa/v0.7.15/bin:$PATH
```
2. Change all *.sh and csmd file to be executable: 
```shell
chmod +x *.sh; chmod +x csmd
```
3. Add csmd into your PATH available: 
```shell
export CSMD_HOME=/public/users/liuyu/csmd_test/code/csmd-master 
export PATH=$PATH:${CSMD_HOME}
```
### 4.3 Databases
CSMD will work with a series of libraries listed in Table 2, including human-related genomes or sequences (21G) and all RefSeq bacteria genomes (150G, as of November 2018). The build process will then require approximately 500GB of additional disk space and 200GB of RAM. These genomes or sequences can be found in DBPATH/hg38/SEQ, DBPATH/AHG/SEQ, DBPATH/EHG/SEQ and DBPATH/RefSeq/bacteria/SEQ, respectively. And the indexed files will be saved in DBPATH/hg38, DBPATH/AHG, DBPATH/EHG and DBPATH/RefSeq/bacteria, respectively.
```shell
csmd --download-library libname --db DBPATH
```
NOTE:<BR/>
--download-library libname: Permissible libname includes “hg38”, “AHG”, “EHG”, “nt”, or “RefSeqBac”.
       --db DBPATH: the store path for the download.
```shell
csmd --build-library libname --db DBPATH
```
NOTE: <BR/>
--build-library: Permissible libname includes “hg38”, “AHG”, “EHG”, “nt”, or “RefSeqBac”. This command generates csmd needed indexed databases, BWA index for hg38 and AHG databases, Bowtie2 index for EHG and RefSeq representative species databases, and blast index for nt and representative genomes databases. To obtain RepBase, go to http://www.girinst.org.<BR/>
--db DBPATH: the store path for the databases, as described above.
### 4.4 Taxonomy files
Taxonomy files include the taxonomic lineage of taxa, information on type strains and material, and host information. This command will download the taxonomy files from NCBI and re-build it according to csmd running. These files can be found in DBPATH/taxonomy/taxdump/ and the re-build file will be saved in DBPATH/taxonomy/taxtree.txt. The files taxdb.tar.gz (ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz) for blast and GenBank2RefSeq.txt in the CSMD download page(/data) for GenBank acceccion number translation should also be saved in DBPATH/taxonomy.
```shell
csmd --download-taxonomy --db DBPATH
csmd --build-taxonomy --db DBPATH
```
NOTE:<BR/>
--download-taxonomy: download taxdump.tar.gz from NCBI.<BR/>
--build-taxonomy: re-organize the taxonomy tree.
## 5. Running CSMD
### 5.1 CS step
This CS procedure performs sensitive and specific computational subtraction of human DNA from the clinical samples. Post-QC paired-end fastq data are assumed as the initial input and indexed human-related genomes or sequences are assumed to be available.
```
csmd --cs hg38Removal --ref REFNAME --thread NUMBER --r1 SAMPLE.R1.fastq --r2 SAMPLE.R2.fastq --output SAMPLE.hg38removal.fastq
csmd --cs ahgRemoval --ref REFNAME --thread NUMBER --input SAMPLE.hg38removal.fastq --output SAMPLE.hg38removal.AHGremoval.fastq
csmd --cs lcrRemoval --ref REFNAME --thread NUMBER --input SAMPLE.hg38removal.AHGremoval.fastq --output SAMPLE.hg38removal.AHGremoval.LCRremoval.fasta
csmd --cs ehgRemoval --ref REFNAME --thread NUMBER \
--input SAMPLE.hg38removal.AHGremoval.LCRremoval.fasta --output SAMPLE.hg38Removal.AHGRemoval.LCRRemoval.EHGRemoval.fasta

```
NOTE:<BR/>
--cs: cs step detail name. Permissible step name includes “hg38Removal”, “ahgRemoval”, “lcrRemoval”, or “ehgRemoval”.<BR/>
--ref: the reference name with the path used for human-derived reads alignment and removal.<BR/>
--thread: number of threads (CPUs) to use.<BR/>
--r1, --r2, --input: the input file with the path. For hg38Removal, the input should be paired-end reads using the parameters “--r1” and “--r2”, but for others, using the parameter “--input". In the phase of hg38Removal, paired-end reads will be combined into single-end data after finishing hg38 reads alignment and removal. The format of the input files is illustrated as the examples.<BR/>
--output: the output file with the path. The format of the output files is illustrated as the examples.<BR/>
### 5.2 MD step
This MD procedure develops a comprehensive and minimally non-redundant reference database using pooled data from the study samples. To overcome the limitations of microbial identification from samples with low microbial biomass and to maximize detection power, all putatively non-human reads in the study group are combined as the input in this step. It includes three key sub-steps to make the microbiome finding as accurate as possible: species finding, species correction and species refinement.
```shell
csmd --md finding --ref REFNAME --thread NUMBER --input pooled_nonHuman.fasta --outdir OUTDIR
csmd --md correction --sam csmd_finding.sam --report csmd_finding_report.tsv \
-–cutoff 25 --thread NUMBER --nt NTNAME --refseq RSNAME --taxdir TAXDIR --outdir OUTDIR
csmd --md refinement --seqdir SEQPATH --seqlist DBLIST --thread NUMBER --input pooled_nonHuman.fasta --outdir OUTDIR
```
### 5.3 Profile step
This procedure provides accurate taxonomy classification for each sample based on a mapping of metagenomic reads against the comprehensive and minimally non-redundant reference database generated in the MD step.
```shell
csmd --pf dbsetup --seqdir SEQPATH --seqlist DBUPDATE --outdir OUTDIR
csmd --pf profile --ref REFNAME --thread NUMBER \
--input SAMPLE.hg38Removal.AHGRemoval.LCRRemoval.EHGRemoval.fasta --sam SAMPLE.csmd.sam –-report SAMPLE.csmd.profile.report
```
NOTE: --pf: profile step detail name. Permissible step name includes “dbsetdup” for the indexed CSMD database and “profile” for single sample microbiome profiling.<BR/>
--seqdir, --seqlist: see MD step.<BR/>
--sam: reads alignment detail in SAM format.<BR/>
--report: pathoscope style report in tsv format.<BR/>
