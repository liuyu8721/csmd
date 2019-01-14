args <- commandArgs(TRUE)
r <- as.character(args[1])
x <- as.integer(args[2])
n <- as.integer(args[3])

calCoverage <- function(file,genomeLength,step) {

  	fileName <- basename(file)

  	if(genomeLength > 0) {
    		if (file.info(file)$size==0) {
    		  count <- rep(0, length(seq(1, genomeLength, step)))
    		} else {
    		  count <- NULL
    		  alldata <- read.csv(file, header = F, stringsAsFactors = FALSE, sep = "\t")
    		  index <- step
    		  while(index <= genomeLength + step) {
      		  count <- c(count, nrow(na.omit(alldata[alldata$V4 > (index - step) & alldata$V4 < index,])))
      			index <- index + step
    		  }
    		}
    		# print(count)
		    write.csv(c(genomeLength, count), paste0(file, ".csv"), row.names = F, quote = F, na = "NA")	
  } else {
    print("ERROR: invalid genome length")
  }
}
  
calCoverage(r, x, n)
