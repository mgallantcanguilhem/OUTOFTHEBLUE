library(qqman)
library(tidyverse)
library(ggbio)
library(GenomicRanges)
library(rtracklayer)  # For importing GTF/GFF files
library(ggpubr)


setwd('~/Documents/CollègeFr/Stage/Metrics/11_plink')

par(mfrow=c(2,2))

# Load data
gwas <- read.table("gwas_results.assoc", header=T)
gwas <- na.omit(gwas)

# Manhattan plot
manhattan(gwas, col=c("goldenrod", "dodgerblue3"), suggestiveline = FALSE, genomewideline = -log10(2.007282e-07), ylim=c(-0.5,8))
mtext("A", side = 3, line = -1.5, adj = 0.05, cex = 1.5, font = 2)

# Load data
gwas <- read.table("gwas_results_brightness.assoc.linear", header=T)
gwas <- na.omit(gwas)

# Manhattan plot
manhattan(gwas, col=c("goldenrod", "dodgerblue3"), suggestiveline = FALSE, genomewideline = -log10(2.007282e-07), ylim=c(-0.5,8))
mtext("B", side = 3, line = -1.5, adj = 0.05, cex = 1.5, font = 2)

# Load data
gwas <- read.table("gwas_results_quanti2_all.assoc.linear", header=T)
gwas <- na.omit(gwas)

# Manhattan plot
manhattan(gwas, col=c("goldenrod", "dodgerblue3"), suggestiveline = FALSE, genomewideline = -log10(2.007282e-07), ylim=c(-0.5,8))
mtext("C", side = 3, line = -1.5, adj = 0.05, cex = 1.5, font = 2)

# Load data
gwas <- read.table("gwas_results_blue_PC1.assoc.linear", header=T)
gwas <- na.omit(gwas)

# Manhattan plot
manhattan(gwas, col=c("goldenrod", "dodgerblue3"), suggestiveline = FALSE, genomewideline = -log10(2.007282e-07), ylim=c(-0.5,8))
mtext("D", side = 3, line = -1.5, adj = 0.05, cex = 1.5, font = 2)

# Load data
gwas <- read.table("gwas_results_quanti_all.assoc.linear", header=T)
gwas <- na.omit(gwas)
gwas <- gwas %>% filter(P<0.05)

# Manhattan plot
manhattan(gwas, col=c("goldenrod", "dodgerblue3"), suggestiveline = FALSE, genomewideline = -log10(2.007282e-07), ylim=c(-0.5,8))

# Load annotation file (GFF/GTF format)
gtf <- rtracklayer::import("Morpho_telemachus_blue.gff3")
# Convert to GenomicRanges object
gr <- GRanges(gtf)

# Filter genes only in Chr30 and the specified range
chr30_genes <- gr[seqnames(gr) == "ptg000030l" & start(gr) >= 500000 & end(gr) <= 2500000]
chr30_genes <- chr30_genes[chr30_genes$type == "gene", ]

gwas_30 <- gwas %>% filter(CHR==30)

# Manhattan plot
p1 <- ggplot(gwas_30, aes(x = BP, y = -log10(P))) +
  geom_point(data = gwas_30[1490000>gwas_30$BP | gwas_30$BP>1550000,], color = "dodgerblue3") +
  geom_point(data = gwas_30[1490000<gwas_30$BP & gwas_30$BP<1550000,], color = "brown3") +  # Green points for specific SNPs
  xlim(1300000,1700000) +
  geom_hline(yintercept = -log10(2.007282e-07), color = "red") +
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12,14), name ="-log10(p)") +
  labs(x = "Genomic Position (Chr30)") +
  theme_classic() +
  theme(axis.text.x = element_blank())

p2 <- autoplot(chr30_genes) + xlim(1250000,1750000) + theme_minimal() + theme(panel.grid = element_blank())
tracks(p1, p2, heights=c(10,1), xlab = "Genomic Position (Chr30)")