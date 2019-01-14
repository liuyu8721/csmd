args <- commandArgs(TRUE)
TAXTREE <- args[1]
WDPATH <- args[2]
THRESHOLD <- as.numeric(args[3])

sname <- read.delim2(paste0(TAXTREE,"/taxtree.txt"), header = T, stringsAsFactors = F)[,1:4]
sname$taxid <- as.double(sname$genome.taxid)
names(sname)[4] <- "species"

setwd(WDPATH)
id.list <- read.delim("screening_species.list", header = F, stringsAsFactors = F)$V1

pool.nt <- NULL
pool.nt.rep <- NULL
for (rsid in id.list) {
  print(rsid)

  #------------------- Results of NT-based alignment -------------------#
  NT <- read.table(paste0("./BLAST/nt/", rsid, ".sample.nt.blast"), header = F, sep = "\t",
                   colClasses = "character",
                   col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                                 "sstart", "send", "evalue", "bitscore", "taxid"))
  for (i in 3:ncol(NT)) {
    NT[,i] <- as.double(NT[,i])
  }

  NT.order <- NT[with(NT, order(qseqid, evalue, - bitscore, - pident, - length)),]

  if (nrow(NT)>0) {
      NT.order.rep <- NULL
      SeqID.list <- unique(NT.order$qseqid)
      for (i in 1:length(SeqID.list)) {
        x <- NT.order[NT.order$qseqid == SeqID.list[i],]
        y <- x[1,]
        z <- x[with(x, evalue==y$evalue, bitscore==y$bitscore, pident==y$pident, length=y$length), ]
        z$weight <- 1 / nrow(z)
        NT.order.rep <- rbind(NT.order.rep, z)
      }

      NT.order.sname <- merge(NT.order, sname, by = "taxid", all.x = T, all.y = F)
      NT.order.rep.sname <- merge(NT.order.rep, sname, by = "taxid", all.x = T, all.y = F)
      NT.order.rep.sname <- NT.order.rep.sname[with(NT.order.rep.sname, order(qseqid)), ]

      ytab <- xtabs(weight ~ species, data = NT.order.rep.sname)
      ytab.det <- data.frame(rsid = rsid,
                             ntotal = length(readLines(paste0("./FASTQ/", rsid, ".fastq"))) / 4,
                             nsampling = length(readLines(paste0("./FASTA/", rsid, ".sample.fasta"))) / 2,
                             nblast = length(SeqID.list),
                             nlink.spe = sum(ytab),
                             spe = names(ytab),
                             spe.freq = as.vector(ytab),
                             stringsAsFactors = FALSE)
      ytab.det <- ytab.det[with(ytab.det, order(- spe.freq)), ]

      ytab.det$str <- NA
      ytab.det$str.freq <- NA
      for (j in 1:length(ytab.det$spe)) {
        x <- NT.order.rep.sname[NT.order.rep.sname$species==ytab.det$spe[j], ]
        str.freq <- xtabs(weight ~ sseqid, data = x)
        ytab.det$str[j] <- names(str.freq)[which(str.freq==max(str.freq))][1]
        ytab.det$str.freq[j] <- max(str.freq)
      }
      print(ytab.det)

      if (rsid==id.list[1]) {
	      pool.nt <- ytab.det
      } else {
        pool.nt <- rbind(pool.nt, ytab.det)
      }
      if (rsid==id.list[1]) {
        if (sum(ytab.det$spe.freq>THRESHOLD)>0) {
		      pool.nt.rep <- ytab.det[which(ytab.det$spe.freq>THRESHOLD),]
	      } else pool.nt.rep <- ytab.det[1,]
      } else {
	      if (sum(ytab.det$spe.freq>THRESHOLD)>0) {
		      pool.nt.rep <- rbind(pool.nt.rep, ytab.det[which(ytab.det$spe.freq>THRESHOLD),])
	      } else pool.nt.rep <- rbind(pool.nt.rep, ytab.det[1,])
      }
  }

  write.table(pool.nt, paste0("./pool.nt_summary.txt"), row.names = F, col.names = T, sep = "\t", quote = F, na = "NA")
  write.table(pool.nt.rep, paste0("./pool.nt.rep_summary.txt"), row.names = F, col.names = T, sep = "\t", quote = F, na = "NA")
}
