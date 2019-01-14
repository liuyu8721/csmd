args <- commandArgs(TRUE)
DBPATH <- args[1]

setwd(DBPATH)
x <- read.delim("assembly_summary.txt", header = F, sep = "\t", skip = 2, stringsAsFactors = F, col.names = c("assembly_accession", "bioproject", "biosample", "wgs_master", "refseq_category", "taxid", "species_taxid", "organism_name", "infraspecific_name", "isolate", "version_status", "assembly_level", "release_type", "genome_rep", "seq_rel_date", "asm_name", "submitter", "gbrs_paired_asm", "paired_asm_comp", "ftp_path", "excluded_from_refseq", "relation_to_type_material"))
x$seq_date <- as.Date(x$seq_rel_date, "%Y/%m/%d")

genomes.taxid <- x$taxid[!duplicated(x$taxid)]
RepGenomes <- NULL
for (i in 1:length(genomes.taxid)) {
  y <- x[which(x$taxid==genomes.taxid[i]),]
  if (nrow(y)==1) {
    RepGenomes <- rbind(RepGenomes, y)
  } else {
    y$L1 <- NA
    y$L1[y$assembly_level=="Complete Genome"] <- 1
    y$L1[y$assembly_level=="Chromosome"] <- 2
    y$L1[y$assembly_level=="Scaffold"] <- 3
    y$L1[y$assembly_level=="Contig"] <- 4
    y$L2 <- NA
    y$L2[y$refseq_category=="reference genome"] <- 1
    y$L2[y$refseq_category=="representative genome"] <- 2
    y$L2[y$refseq_category=="na"] <- 3
    
    y <- y[with(y, order(L1, L2, - as.numeric(seq_date))),]
    RepGenomes <- rbind(RepGenomes, y[1,-which(names(y) %in% c("L1", "L2"))])
  }
}
write.csv(RepGenomes, "GenomeRepresentative.csv", quote = F, row.names = F)
