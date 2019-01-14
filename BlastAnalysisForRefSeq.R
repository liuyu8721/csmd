args <- commandArgs(TRUE)
TAXDIR <- args[1]
WDPATH <- args[2]
THRESHOLD <- as.numeric(args[3])

sname <- read.delim2(paste0(TAXDIR, "/taxtree.txt"), header = T, stringsAsFactors = F)[,1:4]
sname$taxid <- as.double(sname$genome.taxid)
names(sname)[4] <- "species"

GenBank2RefSeq <- read.delim2(paste0(TAXDIR, "/GenBank2RefSeq.txt"), header = T, stringsAsFactors = F)
RefSeq2taxid <- GenBank2RefSeq[, c("taxid", "RefSeq_Accn")]

setwd(WDPATH)
id.list <- read.delim("screening_species.list", header = F, stringsAsFactors = F)$V1

pool.RefSeq <- NULL
pool.RefSeq.rep <- NULL
for (rsid in id.list) {
  print(rsid)
  
  #------------------- Results of RefSeq-based alignment -------------------#
  RefSeq <- read.table(paste0("./BLAST/RefSeq/", rsid, ".sample.RefSeq.blast"), header = F, sep = "\t",
                       colClasses = "character",
                       col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                                     "sstart", "send", "evalue", "bitscore"))
  RefSeq$RefSeq_Accn <- RefSeq$sseqid
  RefSeq <- merge(RefSeq2taxid, RefSeq, by = "RefSeq_Accn", all.x = F, all.y = T, sort = F)
  for (i in 5:ncol(RefSeq)) {
    RefSeq[,i] <- as.double(RefSeq[,i])
  }
  
  RefSeq.order <- RefSeq[with(RefSeq, order(qseqid, evalue, - bitscore, - pident, - length)),]
  
  if (nrow(RefSeq)>0) {
    RefSeq.order.rep <- NULL
    SeqID.list <- unique(RefSeq.order$qseqid)
    for (i in 1:length(SeqID.list)) {
      x <- RefSeq.order[RefSeq.order$qseqid == SeqID.list[i],]
      y <- x[1,]
      z <- x[with(x, evalue==y$evalue, bitscore==y$bitscore, pident==y$pident, length=y$length), ]
      z$weight <- 1 / nrow(z)
      RefSeq.order.rep <- rbind(RefSeq.order.rep, z)
    }
    
    RefSeq.order.sname <- merge(RefSeq.order, sname, by = "taxid", all.x = T, all.y = F)
    RefSeq.order.rep.sname <- merge(RefSeq.order.rep, sname, by = "taxid", all.x = T, all.y = F)
    RefSeq.order.rep.sname <- RefSeq.order.rep.sname[with(RefSeq.order.rep.sname, order(qseqid)), ]
    
    ytab <- xtabs(weight ~ species, data = RefSeq.order.rep.sname)
    ytab.det <- data.frame(rsid = rsid,
                           ntotal = length(readLines(paste0("./FASTQ/", rsid, ".fastq"))) / 4,
                           nsampling = length(readLines(paste0("./FASTA/", rsid, ".sample.fasta"))) / 2,
                           nblast = length(SeqID.list),
                           nlink.spe = sum(ytab),
                           spe = names(ytab),
                           spe.freq = as.vector(ytab),
                           stringsAsFactors = FALSE)
    ytab.det <- ytab.det[with(ytab.det, order(- spe.freq)), ]
    
    ytab.det$str <- rep("", nrow(ytab.det))
    ytab.det$str.freq <- rep(NA, nrow(ytab.det))
    for (j in 1:length(ytab.det$spe)) {
      x <- RefSeq.order.rep.sname[RefSeq.order.rep.sname$species==ytab.det$spe[j], ]
      str.freq <- xtabs(weight ~ sseqid, data = x)
      ytab.det$str[j] <- names(str.freq)[which(str.freq==max(str.freq))][1]
      ytab.det$str.freq[j] <- max(str.freq)
    }
    print(ytab.det)
    
    if (rsid==id.list[1]) {
      pool.RefSeq <- ytab.det
    } else {
      pool.RefSeq <- rbind(pool.RefSeq, ytab.det)
    }
    if (rsid==id.list[1]) {
      if (sum(ytab.det$spe.freq>THRESHOLD)>0) {
        pool.RefSeq.rep <- ytab.det[which(ytab.det$spe.freq>THRESHOLD),]
      } else pool.RefSeq.rep <- ytab.det[1,]
    } else {
      if (sum(ytab.det$spe.freq>THRESHOLD)>0) {
        pool.RefSeq.rep <- rbind(pool.RefSeq.rep, ytab.det[which(ytab.det$spe.freq>THRESHOLD),])
      } else pool.RefSeq.rep <- rbind(pool.RefSeq.rep, ytab.det[1,])
    }
  }
  
  write.table(pool.RefSeq, paste0("./pool.RefSeq_summary.txt"), row.names = F, col.names = T, sep = "\t", quote = F, na = "NA")
  write.table(pool.RefSeq.rep, paste0("./pool.RefSeq.rep_summary.txt"), row.names = F, col.names = T, sep = "\t", quote = F, na = "NA")
}
