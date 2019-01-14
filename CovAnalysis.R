args <- commandArgs(TRUE)
TLPATH <- args[1]
WDPATH <- args[2]

setwd(TLPATH)
#print(getwd())
source("cov.depth.analysis.r")

setwd(WDPATH)
library(stringr)
flist <- list.files()
id.list <- unique(word(flist[grep(".csv", flist)], 1, 2, fixed(".")))

### Bin coverage
cov.stat <- data.frame(id = id.list,
                       BCov = NA,
                       BCov.uniq = NA,
                       mu = NA,
                       pi = NA,
                       sigma = NA,
                       stringsAsFactors = F)

for (i in 1:length(id.list)) {
  
  ### basic info
  id <- id.list[i]
  print(id)

  ### Read data
  x <- read.csv(paste0(id,".csv"), header = T)
  gs <- x[1, 1]
  y <- x[-1, ]
  x.uniq <- read.csv(paste0(id,".unique.csv"), header = T)
  y.uniq <- x.uniq[-1,]

  ### Coverage calculation
  cov.stat$BCov[i] <- (1 - length(which(y==0)) / length(y)) * 100
  cov.stat$BCov.uniq[i] <- (1 - length(which(y.uniq==0)) / length(y.uniq)) * 100
  
  ### Modeling and plotting
  y.uniq[y.uniq < 1] <- NA
  y.uniq[y.uniq > 0] <- 0.1
  if (cov.stat$BCov[i] > 20) {
    if (cov.stat$BCov[i] > 99) m <- "rnor" else m <- "zinb"
    temp <- cov.depth.analysis(y, method = m, q.outliers = runif(1, min=0.99, max=1), plot = FALSE)
    cov.stat[i, c("mu", "pi", "sigma")] <- c(temp$mu, 1 - temp$pi, temp$sigma)
    
    pdf(paste0(id, ".pdf"))
    plot(y, 
         type = "n", cex = 0.6, font.main = 4, main = id,
         xlim = c(0, length(y)), ylim = c(-0.5, max(y) + 1),
         ylab = "Read counts", xlab = "Position of genome", cex.lab = 1,
         xaxt = "n")
    points(1:length(y), y, pch = 16, col = adjustcolor("grey80", alpha.f = 1), cex = 0.6)
    hist.data <- hist(y, breaks = seq(-0.5, max(y) + 0.5, 1), plot = F)
    L <- length(hist.data$breaks)
    rect(rep(0,  L - 1), hist.data$breaks[1:(L-1)], hist.data$density * 800 , hist.data$breaks[2:L],
         border = adjustcolor("gray40", alpha = 0.8), lwd = 1,
         col = adjustcolor("cyan", alpha.f = 0.6))
    par(new = T)
    plot(1:length(y.uniq), y.uniq, type = "h",
         xlim = c(0, length(y)), ylim = c(-0.5, max(y) + 1),
         ylab = "", xlab = "", col = "blue",
         xaxt = "n", yaxt = "n")
    lines(500 * ( (1 - cov.stat$pi[i]) * c(1, rep(0, length(min(y):max(y)) - 1)) + dnbinom(min(y):max(y), mu = cov.stat$mu[i], size = cov.stat$sigma[i]) ), 
          min(y):max(y),
          col = "red", lwd = 2)
    axis(1, at = seq(0, length(y), 100), labels = format(seq(0, length(y), 100) * 5000, scientific = TRUE), cex.axis = 1)
    text(length(y) * 0.6, max(y) * 0.95, labels = substitute(paste("COV"[italic("bin")], " = ", cov), list(cov=paste0(round(cov.stat$BCov[i], digits = 1), "%"))), pos = 4)
    text(length(y) * 0.6, max(y) * 0.95 - max(y) * 0.1, labels = substitute(paste("uCOV"[italic("bin")], " = ", cov), list(cov=paste0(round(cov.stat$BCov.uniq[i], digits = 1), "%"))), pos = 4)
    text(length(y) * 0.6, max(y) * 0.95 - max(y) * 0.2, labels = substitute(paste("GDV"[italic("bin")], " = ", gdv), list(gdv=round(cov.stat$pi[i], 2))), pos = 4)
    text(length(y) * 0.6, max(y) * 0.95 - max(y) * 0.3, labels = substitute(paste("   ", sigma, " "[italic("bin")], "  = ", BDiv), list(BDiv=round(cov.stat$sigma[i], 1))), pos = 4)
    dev.off()
  }

}
write.csv(cov.stat, file = "../CovStat.csv", quote = F, row.names = T)

DB.update <- cov.stat$id[which(cov.stat$BCov.uniq>0.005 & cov.stat$pi > 0.2 & cov.stat$sigma > 0.5)]
write.table(DB.update, "../DB.update", quote = F, col.names = F, row.names = F)
