args <- commandArgs(TRUE)
DBPATH <- args[1]

require(stringr)
setwd(DBPATH)

nodes <- read.delim2("nodes.dmp", header = F, sep = "\t", stringsAsFactors = F)[, c(1, 3, 5)]
names(nodes) <- c("taxid", "f.taxid", "level")
table(nodes$level)

################################ subSpecies or lower level tree ################################
mr0 <- nodes[nodes$level %in% c("subspecies", "varietas", "forma", "no rank"), ]
table(mr0$level)
taxid.list <- data.frame(taxid = unique(mr0$taxid), stringsAsFactors = F)

mr1 <- merge(taxid.list, nodes, by = "taxid", all.x = T, all.y = F, sort = F)
names(mr1) <- c("T0", "F0", "L0")
mr1[mr1$T0==1, "L0"] <- "root"
table(mr1$L0)
mr1 <- mr1[mr1$L0!="root",]

mr2 <- mr1
mr2$taxid <- mr2$F0
mr2 <- merge(mr2, nodes, by = "taxid", all.x = T, all.y = F, sort = F)
mr2$TS <- mr2$F0
mr2$FS <- mr2$f.taxid
mr2$LS <- mr2$level
if (length(which(mr2$TS==1))) mr2[which(mr2$TS==1),"LS"] <- "root"
mr2 <- mr2[, -which(names(mr2) %in% c("taxid", "f.taxid", "level"))]
table(mr2$LS)
mr2.update <- NULL
for (i in 1:17) {
  mr2.update <- rbind(mr2.update, mr2[which(mr2$LS %in% c("species", "genus", "family", "order", "class", "phylum", "superkingdom", "root")),])
  if (length(which(mr2$LS %in% c("species", "genus", "family", "order", "class", "phylum", "superkingdom", "root")))) {
    bad <- mr2[-which(mr2$LS %in% c("species", "genus", "family", "order", "class", "phylum", "superkingdom", "root")),]
  } else bad <- mr2
  bad$taxid <- bad$FS
  bad <- merge(bad, nodes, by = "taxid", all.x = T, all.y = F, sort = F)[,-1]
  bad$TS <- bad$FS
  bad$FS <- bad$f.taxid
  bad$LS <- bad$level
  # if (length(which(bad$TG==1))) bad[which(bad$TG==1),"LG"] <- "root"
  mr2 <- bad[,-which(names(bad) %in% c("f.taxid", "level"))]
  if (length(which(mr2$TS==1))) mr2[which(mr2$TS==1),"LS"] <- "root"
  print(table(mr2$LS))
}
table(mr2.update$LS)

# mr.1 <- mr2.update[which(mr2.update$LS=="species"),]
# mr.2 <- mr2.update[-which(mr2.update$LS=="species"),]

################################ Species tree ################################
m0 <- nodes[nodes$level %in% c("species"), ]
table(m0$level)
taxid.list <- data.frame(taxid = unique(m0$taxid), stringsAsFactors = F)

m1 <- merge(taxid.list, nodes, by = "taxid", all.x = T, all.y = F, sort = F)
names(m1) <- c("T0", "F0", "L0")
m1$TS <- m1$T0
m1$FS <- m1$F0
m1$LS <- m1$L0
m1.update <- rbind(mr2.update, m1)
table(m1.update$LS)

m2.end <- m1.update[which(m1.update$LS %in% c("root")),]
m2 <- m1.update[-which(m1.update$LS %in% c("root")),]
m2$taxid <- m2$FS
m2 <- merge(m2, nodes, by = "taxid", all.x = T, all.y = F, sort = F)
m2$TG <- m2$FS
m2$FG <- m2$f.taxid
m2$LG <- m2$level
if (length(which(m2$TG==1))) m2[which(m2$TG==1),"LG"] <- "root"
m2 <- m2[, -which(names(m2) %in% c("taxid", "f.taxid", "level"))]
table(m2$LG)
m2.update <- NULL
for (i in 1:18) {
  m2.update <- rbind(m2.update, m2[which(m2$LG %in% c("genus", "family", "order", "class", "phylum", "superkingdom", "root")),])
  if (length(which(m2$LG %in% c("genus", "family", "order", "class", "phylum", "superkingdom", "root")))) {
    bad <- m2[-which(m2$LG %in% c("genus", "family", "order", "class", "phylum", "superkingdom", "root")),]
  } else bad <- m2
  bad$taxid <- bad$FG
  bad <- merge(bad, nodes, by = "taxid", all.x = T, all.y = F, sort = F)[,-1]
  bad$TG <- bad$FG
  bad$FG <- bad$f.taxid
  bad$LG <- bad$level
  # if (length(which(bad$TG==1))) bad[which(bad$TG==1),"LG"] <- "root"
  m2 <- bad[,-which(names(bad) %in% c("f.taxid", "level"))]
  if (length(which(m2$TG==1))) m2[which(m2$TG==1),"LG"] <- "root"
  print(table(m2$LG))
}
table(m2.update$LG)

m3.end <- m2.update[which(m2.update$LG %in% c("root")),]
m3 <- m2.update[-which(m2.update$LG %in% c("root")),]
m3$taxid <- m3$FG
m3 <- merge(m3, nodes, by = "taxid", all.x = T, all.y = F, sort = F)
m3$TF <- m3$FG
m3$FF <- m3$f.taxid
m3$LF <- m3$level
if (length(which(m3$TF==1))) m3[which(m3$TF==1),"LF"] <- "root"
m3 <- m3[, -which(names(m3) %in% c("taxid", "f.taxid", "level"))]
table(m3$LF)
m3.update <- NULL
for (i in 1:18) {
  m3.update <- rbind(m3.update, m3[which(m3$LF %in% c("family", "order", "class", "phylum", "superkingdom", "root")),])
  if (length(which(m3$LF %in% c("family", "order", "class", "phylum", "superkingdom", "root")))) {
    bad <- m3[-which(m3$LF %in% c("family", "order", "class", "phylum", "superkingdom", "root")),]
  } else bad <- m3
  bad$taxid <- bad$FF
  bad <- merge(bad, nodes, by = "taxid", all.x = T, all.y = F, sort = F)[,-1]
  bad$TF <- bad$FF
  bad$FF <- bad$f.taxid
  bad$LF <- bad$level
  m3 <- bad[,-which(names(bad) %in% c("f.taxid", "level"))]
  if (length(which(m3$TF==1))) m3[which(m3$TF==1),"LF"] <- "root"
  print(table(m3$LF))
}
table(m3.update$LF)

m4.end <- m3.update[which(m3.update$LF %in% c("root")),]
m4 <- m3.update[-which(m3.update$LF %in% c("root")),]
m4$taxid <- m4$FF
m4 <- merge(m4, nodes, by = "taxid", all.x = T, all.y = F, sort = F)
m4$TO <- m4$FF
m4$FO <- m4$f.taxid
m4$LO <- m4$level
if (length(which(m4$TO==1))) m4[which(m4$TO==1),"LO"] <- "root"
m4 <- m4[, -which(names(m4) %in% c("taxid", "f.taxid", "level"))]
table(m4$LO)
m4.update <- NULL
for (i in 1:18) {
  m4.update <- rbind(m4.update, m4[which(m4$LO %in% c("order", "class", "phylum", "superkingdom", "root")),])
  if (length(which(m4$LO %in% c("order", "class", "phylum", "superkingdom", "root")))) {
    bad <- m4[-which(m4$LO %in% c("order", "class", "phylum", "superkingdom", "root")),]
  } else bad <- m4
  bad$taxid <- bad$FO
  bad <- merge(bad, nodes, by = "taxid", all.x = T, all.y = F, sort = F)[,-1]
  bad$TO <- bad$FO
  bad$FO <- bad$f.taxid
  bad$LO <- bad$level
  m4 <- bad[,-which(names(bad) %in% c("f.taxid", "level"))]
  if (length(which(m4$TO==1))) m4[which(m4$TO==1),"LO"] <- "root"
  print(table(m4$LO))
}
table(m4.update$LO)

m5.end <- m4.update[which(m4.update$LO %in% c("root")),]
m5 <- m4.update[-which(m4.update$LO %in% c("root")),]
m5$taxid <- m5$FO
m5 <- merge(m5, nodes, by = "taxid", all.x = T, all.y = F, sort = F)
m5$TC <- m5$FO
m5$FC <- m5$f.taxid
m5$LC <- m5$level
if (length(which(m5$TC==1))) m5[which(m5$TC==1),"LC"] <- "root"
m5 <- m5[, -which(names(m5) %in% c("taxid", "f.taxid", "level"))]
table(m5$LC)
m5.update <- NULL
for (i in 1:18) {
  m5.update <- rbind(m5.update, m5[which(m5$LC %in% c("class", "phylum", "superkingdom", "root")),])
  if (length(which(m5$LC %in% c("class", "phylum", "superkingdom", "root")))) {
    bad <- m5[-which(m5$LC %in% c("class", "phylum", "superkingdom", "root")),]
  } else bad <- m5
  bad$taxid <- bad$FC
  bad <- merge(bad, nodes, by = "taxid", all.x = T, all.y = F, sort = F)[,-1]
  bad$TC <- bad$FC
  bad$FC <- bad$f.taxid
  bad$LC <- bad$level
  m5 <- bad[,-which(names(bad) %in% c("f.taxid", "level"))]
  if (length(which(m5$TC==1))) m5[which(m5$TC==1),"LC"] <- "root"
  print(table(m5$LC))
}
table(m5.update$LC)

m6.end <- m5.update[which(m5.update$LC %in% c("root")),]
m6 <- m5.update[-which(m5.update$LC %in% c("root")),]
m6$taxid <- m6$FC
m6 <- merge(m6, nodes, by = "taxid", all.x = T, all.y = F, sort = F)
m6$TP <- m6$FC
m6$FP <- m6$f.taxid
m6$LP <- m6$level
if (length(which(m6$TP==1))) m6[which(m6$TP==1),"LP"] <- "root"
m6 <- m6[, -which(names(m6) %in% c("taxid", "f.taxid", "level"))]
table(m6$LP)
m6.update <- NULL
for (i in 1:18) {
  m6.update <- rbind(m6.update, m6[which(m6$LP %in% c("phylum", "superkingdom", "root")),])
  if (length(which(m6$LP %in% c("phylum", "superkingdom", "root")))) {
    bad <- m6[-which(m6$LP %in% c("phylum", "superkingdom", "root")),]
  } else bad <- m6
  bad$taxid <- bad$FP
  bad <- merge(bad, nodes, by = "taxid", all.x = T, all.y = F, sort = F)[,-1]
  bad$TP <- bad$FP
  bad$FP <- bad$f.taxid
  bad$LP <- bad$level
  m6 <- bad[,-which(names(bad) %in% c("f.taxid", "level"))]
  if (length(which(m6$TP==1))) m6[which(m6$TP==1),"LP"] <- "root"
  print(table(m6$LP))
}
table(m6.update$LP)

m7.end <- m6.update[which(m6.update$LP %in% c("root")),]
m7 <- m6.update[-which(m6.update$LP %in% c("root")),]
m7$taxid <- m7$FP
m7 <- merge(m7, nodes, by = "taxid", all.x = T, all.y = F, sort = F)
m7$TK <- m7$FP
m7$FK <- m7$f.taxid
m7$LK <- m7$level
if (length(which(m7$TK==1))) m7[which(m7$TK==1),"LK"] <- "root"
m7 <- m7[, -which(names(m7) %in% c("taxid", "f.taxid", "level"))]
table(m7$LK)
m7.update <- NULL
for (i in 1:8) {
  m7.update <- rbind(m7.update, m7[which(m7$LK %in% c("superkingdom", "root")),])
  if (length(which(m7$LK %in% c("superkingdom", "root")))) {
    bad <- m7[-which(m7$LK %in% c("superkingdom", "root")),]
  } else bad <- m7
  bad$taxid <- bad$FK
  bad <- merge(bad, nodes, by = "taxid", all.x = T, all.y = F, sort = F)[,-1]
  bad$TK <- bad$FK
  bad$FK <- bad$f.taxid
  bad$LK <- bad$level
  m7 <- bad[,-which(names(bad) %in% c("f.taxid", "level"))]
  if (length(which(m7$TK==1))) m7[which(m7$TK==1),"LK"] <- "root"
  print(table(m7$LK))
}
table(m7.update$LK)

m8.end <- m7.update[which(m7.update$LK %in% c("root")),]
if (length(which(m7.update$LK %in% c("root")))) {
  m8 <- m7.update[-which(m7.update$LK %in% c("root")),]
} else m8 <- m7.update
m8$taxid <- m8$FK
m8 <- merge(m8, nodes, by = "taxid", all.x = T, all.y = F, sort = F)
m8$TR <- m8$FK
m8$FR <- m8$f.taxid
m8$LR <- m8$level
if (length(which(m8$TR==1))) m8[which(m8$TR==1),"LR"] <- "root"
m8 <- m8[, -which(names(m8) %in% c("taxid", "f.taxid", "level"))]
table(m8$LR)
m8.update <- NULL
for (i in 1:2) {
  m8.update <- rbind(m8.update, m8[which(m8$LR %in% c("root")),])
  if (length(which(m8$LR %in% c("root")))) {
    bad <- m8[-which(m8$LR %in% c("root")),]
  } else bad <- m8
  bad$taxid <- bad$FR
  bad <- merge(bad, nodes, by = "taxid", all.x = T, all.y = F, sort = F)[,-1]
  bad$TR <- bad$FR
  bad$FR <- bad$f.taxid
  bad$LR <- bad$level
  m8 <- bad[,-which(names(bad) %in% c("f.taxid", "level"))]
  if (length(which(m8$TR==1))) m8[which(m8$TR==1),"LR"] <- "root"
  print(table(m8$LR))
}
table(m8.update$LR)

all <- merge(m2.end, m3.end, all.x = T, all.y = T, sort = F)
all <- merge(all, m4.end, all.x = T, all.y = T, sort = F)
all <- merge(all, m5.end, all.x = T, all.y = T, sort = F)
all <- merge(all, m6.end, all.x = T, all.y = T, sort = F)
all <- merge(all, m7.end, all.x = T, all.y = T, sort = F)
all <- merge(all, m8.end, all.x = T, all.y = T, sort = F)
all <- merge(all, m8.update, all.x = T, all.y = T, sort = F)
names(all)

# all.update <- all
# for (i in 1:nrow(all.update)) {
#   if (is.na(all.update$LR[i])) {
#     loc <- which(all.update[i, c("LS", "LG", "LF", "LO", "LC", "LP", "LK")] %in% "root")
#     all.update[i, c("TR", "FR", "LR")] <- all.update[i, 1:3 + loc * 3]
#     all.update[i, 1:3 + loc * 3] <- c(NA, NA, NA)
#   }
#   if (is.na(all.update$LK[i])) {
#     loc <- which(all.update[i, c("LS", "LG", "LF", "LO", "LC", "LP")] %in% "superkingdom")
#     if (length(loc)) {
#       all.update[i, c("TK", "FK", "LK")] <- all.update[i, 1:3 + loc * 3]
#       all.update[i, 1:3 + loc * 3] <- c(NA, NA, NA)
#     }
#   }
#   if (is.na(all.update$LP[i])) {
#     loc <- which(all.update[i, c("LS", "LG", "LF", "LO", "LC")] %in% "phylum")
#     if (length(loc)) {
#       all.update[i, c("TP", "FP", "LP")] <- all.update[i, 1:3 + loc * 3]
#       all.update[i, 1:3 + loc * 3] <- c(NA, NA, NA)      
#     }
#   }
#   if (is.na(all.update$LC[i])) {
#     loc <- which(all.update[i, c("LS", "LG", "LF", "LO")] %in% "class")
#     if (length(loc)) {
#       all.update[i, c("TC", "FC", "LC")] <- all.update[i, 1:3 + loc * 3]
#       all.update[i, 1:3 + loc * 3] <- c(NA, NA, NA)      
#     }
#   }
#   if (is.na(all.update$LO[i])) {
#     loc <- which(all.update[i, c("LS", "LG", "LF")] %in% "order")
#     if (length(loc)) {
#       all.update[i, c("TO", "FO", "LO")] <- all.update[i, 1:3 + loc * 3]
#       all.update[i, 1:3 + loc * 3] <- c(NA, NA, NA)      
#     }
#   }
#   if (is.na(all.update$LF[i])) {
#     loc <- which(all.update[i, c("LS", "LG")] %in% "family")
#     if (length(loc)) {
#       all.update[i, c("TF", "FF", "LF")] <- all.update[i, 1:3 + loc * 3]
#       all.update[i, 1:3 + loc * 3] <- c(NA, NA, NA)      
#     }
#   }
#   if (is.na(all.update$LG[i])) {
#     loc <- which(all.update[i, c("LS")] %in% "genus")
#     if (length(loc)) {
#       all.update[i, c("TG", "FG", "LG")] <- all.update[i, 1:3 + loc * 3]
#       all.update[i, 1:3 + loc * 3] <- c(NA, NA, NA)      
#     }
#   }
# }

seq.list <- seq(1, nrow(all), 10000)
for (j in 1:length(seq.list)) {
  all.update <- all[seq.list[j]:min(seq.list[j] + 10000 - 1, nrow(all)),]
  for (i in 1:nrow(all.update)) {
    if (is.na(all.update$LR[i])) {
      loc <- which(all.update[i, c("LS", "LG", "LF", "LO", "LC", "LP", "LK")] %in% "root")
      all.update[i, c("TR", "FR", "LR")] <- all.update[i, 1:3 + loc * 3]
      all.update[i, 1:3 + loc * 3] <- c(NA, NA, NA)
    }
    if (is.na(all.update$LK[i])) {
      loc <- which(all.update[i, c("LS", "LG", "LF", "LO", "LC", "LP")] %in% "superkingdom")
      if (length(loc)) {
        all.update[i, c("TK", "FK", "LK")] <- all.update[i, 1:3 + loc * 3]
        all.update[i, 1:3 + loc * 3] <- c(NA, NA, NA)
      }
    }
    if (is.na(all.update$LP[i])) {
      loc <- which(all.update[i, c("LS", "LG", "LF", "LO", "LC")] %in% "phylum")
      if (length(loc)) {
        all.update[i, c("TP", "FP", "LP")] <- all.update[i, 1:3 + loc * 3]
        all.update[i, 1:3 + loc * 3] <- c(NA, NA, NA)      
      }
    }
    if (is.na(all.update$LC[i])) {
      loc <- which(all.update[i, c("LS", "LG", "LF", "LO")] %in% "class")
      if (length(loc)) {
        all.update[i, c("TC", "FC", "LC")] <- all.update[i, 1:3 + loc * 3]
        all.update[i, 1:3 + loc * 3] <- c(NA, NA, NA)      
      }
    }
    if (is.na(all.update$LO[i])) {
      loc <- which(all.update[i, c("LS", "LG", "LF")] %in% "order")
      if (length(loc)) {
        all.update[i, c("TO", "FO", "LO")] <- all.update[i, 1:3 + loc * 3]
        all.update[i, 1:3 + loc * 3] <- c(NA, NA, NA)      
      }
    }
    if (is.na(all.update$LF[i])) {
      loc <- which(all.update[i, c("LS", "LG")] %in% "family")
      if (length(loc)) {
        all.update[i, c("TF", "FF", "LF")] <- all.update[i, 1:3 + loc * 3]
        all.update[i, 1:3 + loc * 3] <- c(NA, NA, NA)      
      }
    }
    if (is.na(all.update$LG[i])) {
      loc <- which(all.update[i, c("LS")] %in% "genus")
      if (length(loc)) {
        all.update[i, c("TG", "FG", "LG")] <- all.update[i, 1:3 + loc * 3]
        all.update[i, 1:3 + loc * 3] <- c(NA, NA, NA)      
      }
    }
  }
  write.table(all.update, paste0("tree", j, ".txt"), row.names = F, col.names = T, sep = "\t", quote = F, na = "NA")
}

flist <- list.files()
flist <- flist[grep("tree", flist, fixed=TRUE)]
all.update <- NULL
for (i in 1:length(flist)) {
  x <- read.delim2(flist[i], header = T, stringsAsFactors = F)
  all.update <- rbind(all.update, x)
}
file.remove(flist)

snames <- read.delim2("names.dmp", header = F, sep = "\t", quote = "", stringsAsFactors = F)[, seq(1,7,2)]
snames <- snames[which(snames$V7=="scientific name"),][,1:2]
names(snames) <- c("taxid", "sname")
class(snames$taxid)

total <- data.frame(taxid = all.update$T0,
                    T0 = all.update$T0,
                    genome.taxid = all.update$T0,
                    stringsAsFactors = F)
total <- merge(total, snames, by = "taxid", all.x = T, all.y = F, sort = F)
names(total)[4] <- "genome.sname"
total <- merge(total, all.update[,c("T0", "TS")], all.x = T, all.y = T, sort = F)
total$taxid <- total$TS
total <- merge(total, snames, by = "taxid", all.x = T, all.y = F, sort = F)
names(total)[5:6] <- c("species.taxid", "species.sname")
total <- merge(total, all.update[,c("T0", "TG")], all.x = T, all.y = T, sort = F)
total$taxid <- total$TG
total <- merge(total, snames, by = "taxid", all.x = T, all.y = F, sort = F)
names(total)[7:8] <- c("genus.taxid", "genus.sname")
total <- merge(total, all.update[,c("T0", "TF")], all.x = T, all.y = T, sort = F)
total$taxid <- total$TF
total <- merge(total, snames, by = "taxid", all.x = T, all.y = F, sort = F)
names(total)[9:10] <- c("family.taxid", "family.sname")
total <- merge(total, all.update[,c("T0", "TO")], all.x = T, all.y = T, sort = F)
total$taxid <- total$TO
total <- merge(total, snames, by = "taxid", all.x = T, all.y = F, sort = F)
names(total)[11:12] <- c("order.taxid", "order.sname")
total <- merge(total, all.update[,c("T0", "TC")], all.x = T, all.y = T, sort = F)
total$taxid <- total$TC
total <- merge(total, snames, by = "taxid", all.x = T, all.y = F, sort = F)
names(total)[13:14] <- c("class.taxid", "class.sname")
total <- merge(total, all.update[,c("T0", "TP")], all.x = T, all.y = T, sort = F)
total$taxid <- total$TP
total <- merge(total, snames, by = "taxid", all.x = T, all.y = F, sort = F)
names(total)[15:16] <- c("phylum.taxid", "phylum.sname")
total <- merge(total, all.update[,c("T0", "TK")], all.x = T, all.y = T, sort = F)
total$taxid <- total$TK
total <- merge(total, snames, by = "taxid", all.x = T, all.y = F, sort = F)
names(total)[17:18] <- c("superkingdom.taxid", "superkingdom.sname")
total <- merge(total, all.update[,c("T0", "TR")], all.x = T, all.y = T, sort = F)
total$taxid <- total$TR
total <- merge(total, snames, by = "taxid", all.x = T, all.y = F, sort = F)
names(total)[19:20] <- c("root.taxid", "root.sname")
table(total$superkingdom.sname)

write.table(total[,-c(1:2)], "../taxtree.txt", row.names = F, col.names = T, sep = "\t", quote = F, na = "NA")
