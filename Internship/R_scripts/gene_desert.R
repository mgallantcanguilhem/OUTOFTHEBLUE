library(rtracklayer)  # For importing GTF/GFF files
library(tidyverse)
library(ggplot2)

setwd('~/Documents/CollègeFr/Stage/Metrics/gene_desert')

# Load annotation file (GFF/GTF format)
gtf <- rtracklayer::import("Morpho_telemachus_blue.gff3")

# Filter genes only 
gtf_genes <- as.tibble(gtf[gtf$type == "gene", ])

df <- tibble(CHR = 1, POS = gtf_genes$end[1]+1, Width = gtf_genes$start[2] - gtf_genes$end[1])
k=1
for (i in 2:(nrow(gtf_genes)-1)){
  if (gtf_genes$seqnames[i] == gtf_genes$seqnames[i+1]){
    row_supp <- c(paste(k), gtf_genes$end[i]+1, gtf_genes$start[i+1] - gtf_genes$end[i])
    df <- rbind(df, row_supp)
  } else {
    k <- k+1
  }
}

df$CHR <- as.integer(df$CHR)
df$POS <- as.integer(df$POS)
df$Width <- as.integer(df$Width)
df <- df[df$Width>0,]
summary(df)

quantile(df$Width, probs=seq(0.9, 1, 0.01))
df %>% filter(CHR ==30, POS>1400000, POS<1600000)

ggplot(df, aes(x = Width/1000)) + 
  geom_density(fill="dodgerblue2", alpha=0.8) + 
  geom_vline(aes(xintercept=135872/1000), color="red", linetype="dashed", size=0.5) + 
  geom_vline(aes(xintercept=129609/1000), color="blue", linetype="dashed", size=1) +
  xlab("Size of intergenic regions (kb)") + 
  ylab("Density") +
  xlim(c(0,200)) +
  theme_classic()

ggplot(df, aes(x = POS)) + geom_density()

desert <- aggregate(df$Width, by=list(CHR=df$CHR), FUN=sum)
desert$Size <- rep(0, nrow(desert))
k=1
for (i in 1:(nrow(df)-1)){
  if (df$CHR[i]!=df$CHR[i+1]){
    desert$Size[k] <- df$POS[i] + df$Width
    k <- k+1
  }
}
desert$CHR <- as.factor(desert$CHR)
desert$fraction <- 100*desert$x/desert$Size

row.names.remove <- c("27","31","33","34","35","36","37","38","39","40","41")
desert[!(row.names(desert) %in% row.names.remove), ]

ggplot(desert[!(row.names(desert) %in% row.names.remove), ], aes(x=CHR, y=Size, fill=fraction)) + geom_col()
