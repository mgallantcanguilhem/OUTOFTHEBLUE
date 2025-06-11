library(tidyverse)
library(ggpubr)

setwd("~/Documents/CollègeFr/Stage/Metrics/09ter_filter")

var_qual <- read_delim("./mortel_subset.lqual", delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1)
var_qual <- na.omit(var_qual)
var_qual <- var_qual[which(var_qual$qual!=0),]
var_qual <- var_qual[var_qual$qual < 3000, ]

a <- ggplot(var_qual, aes(qual)) +
  geom_density(fill = "dodgerblue2", colour = "black", alpha = 0.3) +
  geom_vline(xintercept = 25, linetype = "dashed", color = "brown2", linewidth=1.2) +
  ylab("Density") +
  xlab("Quality (Phred)") +
  xlim(0,300) +
  theme_classic()

var_depth <- read_delim("./mortel_subset.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
var_depth<- na.omit(var_depth)
var_depth <- var_depth[which(var_depth$mean_depth!=0),]
var_depth <- var_depth[var_depth$mean_depth < 100, ]

b <- ggplot(var_depth, aes(mean_depth)) +
  geom_density(fill = "dodgerblue2", colour = "black", alpha = 0.3) +
  geom_vline(xintercept = c(7,40), linetype = "dashed", color = "brown2", linewidth=1.2) +
  ylab("Density") +
  xlab("Mean depth") +
  xlim(0,70) +
  theme_classic()

var_miss <- read_delim("./mortel_subset.lmiss", delim = "\t",
                      col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

var_miss <- var_miss[which(var_miss$fmiss!=1),]

c <- ggplot(var_miss, aes(fmiss)) +
  geom_density(fill = "dodgerblue2", colour = "black", alpha = 0.3) +
  geom_vline(xintercept = 0.15, linetype = "dashed", color = "brown2", linewidth=1.2) +
  ylab("Density") +
  xlab("Missing frequency") +
  xlim(0,1) +
  theme_classic()

var_freq <- read_delim("./mortel_subset.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 164)
var_freq <- na.omit(var_freq)
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))

d <- ggplot(var_freq, aes(maf)) +
  geom_density(fill = "dodgerblue2", colour = "black", alpha = 0.3) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "brown2", linewidth=1.2) +
  ylab("Density") +
  xlab("Minor allele frequency") +
  theme_classic()

ggarrange(a, b, c, d, labels=c("A", "B", "C", "D"), align = "hv")


ind_depth <- read_delim("./mortel_subset.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)
e <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
e + theme_light()

ind_miss  <- read_delim("./mortel_subset.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
f <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
f + theme_light()

ind_het <- read_delim("./mortel_subset.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
g <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
g + theme_light()

setwd("~/Documents/CollègeFr/Stage/Metrics/09ter_filter/old_gatk")

var_qual <- read_delim("./mortel_subset_old.lqual", delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + geom_vline(xintercept=30)
a + theme_light() + xlim(1,250)

var_depth <- read_delim("./mortel_subset_old.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
b <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + geom_vline(xintercept=c(7,40))
b + theme_light() + xlim(1,70)

var_miss <- read_delim("./mortel_subset_old.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
c <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + geom_vline(xintercept=0.15)
c + theme_light() + xlim(0.01, 0.9)

var_freq <- read_delim("./mortel_subset_old.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 164)
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
d <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + geom_vline(xintercept=0.1)
d + theme_light() + xlim(0, 0.5)

ind_depth <- read_delim("./mortel_subset_old.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)
e <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
e + theme_light()

ind_miss  <- read_delim("./mortel_subset_old.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
f <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
f + theme_light()

ind_het <- read_delim("./mortel_subset_old.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
g <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
g + theme_light()
