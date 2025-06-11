#### coverage profile by individual ####

#set to pass arguments from bash
args = commandArgs(trailingOnly=TRUE)

dpt <- read.table(args[1])

#add colnames
colnames(dpt) <- c("chr","pos", paste("ind",c(1:(ncol(dpt)-2)), sep=""))

#create a new position genome wide
POS_new=NULL
position_max = 0
for (i in levels(as.factor(dpt$chr))){
  pos_relative = subset(dpt, chr == i)$pos
  POS_new = c(POS_new, pos_relative + position_max)
  position_max = max(pos_relative)+position_max
}
dpt$cum_pos = POS_new

#new df
dptind <- as.data.frame(dpt$cum_pos)
colnames(dptind) <- "cum_pos"

# Compute bin number for sliding window
dptind$bin <- (dptind$cum_pos - 1) %/% 100000 + 1
depthdf <- aggregate(cum_pos ~ bin, data = dptind, FUN = mean)
colnames(depthdf)[2] <- "midpos"
depthdf$midpos <- round(depthdf$midpos,digits = 0)

for (i in 1:(ncol(dpt)-2)) {
  #select individual
  ind <- dpt[,(2+i)]
  # Compute mean for each bin
  means <- as.data.frame(tapply(ind, dptind$bin, mean))
  #change colname according to original col
  colnames(means) <- colnames(dpt)[i+2]
  #merge with general df
  depthdf <- cbind(depthdf,means)
}

write.table(depthdf, file = args[2], quote=FALSE, row.names=FALSE)
