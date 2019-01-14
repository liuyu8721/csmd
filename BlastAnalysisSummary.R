args <- commandArgs(TRUE)
TAXDIR <- args[1]
WDPATH <- args[2]
THRESHOLD <- as.numeric(args[3])

GenBank2RefSeq <- read.delim2(paste0(TAXDIR, "/GenBank2RefSeq.txt"), header = T, stringsAsFactors = F)
RefSeq2taxid <- GenBank2RefSeq[, c("taxid", "RefSeq_Accn")]
GenBank2RefSeq.uniq <- GenBank2RefSeq[which(GenBank2RefSeq$Flag==1), ]
GenBank2RefSeq.uniq <- GenBank2RefSeq.uniq[-which(duplicated(GenBank2RefSeq$RefSeqUAccn)),c("RefSeqUAccn", "species.sname")]
names(GenBank2RefSeq.uniq)[1] <- "RefSeq_Accn.species"

setwd(WDPATH)
pool.nt.rep <- read.delim2("pool.nt.rep_summary.txt", header = T, stringsAsFactors = F)
pool.nt.rep$spe.freq <- as.numeric(pool.nt.rep$spe.freq)
pool.RefSeq.rep <- read.delim2("pool.RefSeq.rep_summary.txt", header = T, stringsAsFactors = F)
pool.RefSeq.rep$spe.freq <- as.numeric(pool.RefSeq.rep$spe.freq)
pool.nt.rep.update <- NULL
for (i in 1:nrow(pool.nt.rep)) {
 id <- pool.nt.rep$rsid[i]
 if (pool.nt.rep$spe.freq[i] <= THRESHOLD) {
   temp <- pool.RefSeq.rep[which(pool.RefSeq.rep$rsid == id & pool.RefSeq.rep$spe.freq > THRESHOLD),]
   if (length(nrow(temp)) > 0) pool.nt.rep.update <- rbind(pool.nt.rep.update, temp)
 } else pool.nt.rep.update <- rbind(pool.nt.rep.update, pool.nt.rep[i, ])
}

pool.nt.rep.update$GenBank_Accn <- pool.nt.rep.update$str
pool.nt.rep.update <- merge(pool.nt.rep.update, GenBank2RefSeq[,c("GenBank_Accn", "RefSeqUAccn")], by = "GenBank_Accn", all.x = T, all.y = F)
pool.nt.rep.update$RefSeq_Accn <- pool.nt.rep.update$str
pool.nt.rep.update <- merge(pool.nt.rep.update, GenBank2RefSeq[,c("RefSeq_Accn", "RefSeqUAccn")], by = "RefSeq_Accn", all.x = T, all.y = F)
pool.nt.rep.update$species.sname <- pool.nt.rep.update$spe
pool.nt.rep.update <- merge(pool.nt.rep.update, GenBank2RefSeq.uniq, by = "species.sname", all.x = T, all.y = F)

pool.nt.rep.update$RefSeq.AA <- NA
pool.nt.rep.update$RefSeq.AA[which(!is.na(pool.nt.rep.update$RefSeqUAccn.x))] <- pool.nt.rep.update$RefSeqUAccn.x[which(!is.na(pool.nt.rep.update$RefSeqUAccn.x))]
pool.nt.rep.update$RefSeq.AA[which(!is.na(pool.nt.rep.update$RefSeqUAccn.y))] <- pool.nt.rep.update$RefSeqUAccn.y[which(!is.na(pool.nt.rep.update$RefSeqUAccn.y))]
pool.nt.rep.update$RefSeq.AA[which(is.na(pool.nt.rep.update$RefSeq.AA))] <- pool.nt.rep.update$RefSeq_Accn.species[which(is.na(pool.nt.rep.update$RefSeq.AA))]
write.table(pool.nt.rep.update[, c("rsid", "ntotal", "nsampling", "nblast", "nlink.spe", "spe", "spe.freq", "str", "str.freq", "RefSeqUAccn.x", "RefSeqUAccn.y", "RefSeq_Accn.species", "RefSeq.AA")], "pool.nt.rep_update.txt", row.names = F, col.names = T, sep = "\t", quote = F, na = "NA")

pool.rep.update <- pool.nt.rep.update[which(!is.na(pool.nt.rep.update$RefSeq.AA)), !(names(pool.nt.rep.update) %in% c("RefSeq_Accn", "GenBank_Accn", "RefSeqUAccn.x", "RefSeqUAccn.y", "RefSeq_Accn.species", "species.sname"))]
write.table(pool.rep.update, "pool.rep_update.txt", row.names = F, col.names = T, sep = "\t", quote = F, na = "NA")

pool.rep.update$str.freq <- as.numeric(pool.rep.update$str.freq)
DB.list <- pool.rep.update[with(pool.rep.update, order(spe, - spe.freq, - str.freq)),]
DB.list <- DB.list[!duplicated(DB.list$spe), "RefSeq.AA"]
write.table(DB.list, "DB.list", row.names = F, col.names = F, sep = "\t", quote = F, na = "NA")

