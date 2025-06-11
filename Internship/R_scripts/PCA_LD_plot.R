library(tidyverse)
library(ggplot2)
library(rstatix)

setwd('~/Documents/CollègeFr/Stage/Metrics/11_plink')

# read in data
pca <- read_table2("mortel.eigenvec", col_names = FALSE)
eigenval <- scan("mortel.eigenval")

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# sort out the individual location and color
# color
col <- rep(NA, length(pca$ind))
col[grep("blue", pca$ind)] <- "Blue"
col[grep("oran", pca$ind)] <- "Orange"
# location
loc <- rep(NA, length(pca$ind))
loc[grep("CAC", pca$ind)] <- "Cacao"
loc[grep("KAW", pca$ind)] <- "Kaw-Patawa"
loc[grep("PAT", pca$ind)] <- "Kaw-Patawa"
# combine - if you want to plot each in different colours
col_loc <- paste0(col, "_", loc)

# remake data.frame
pca <- as_tibble(data.frame(pca, col, loc, col_loc))

# first convert to percentage variance explained
pve <- data.frame(PC = 1:length(eigenval), pve = eigenval/sum(eigenval)*100)
# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_classic()

# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = col, shape = loc)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("blue", "orange"), name="Color") + scale_shape(name="Location")
b <- b + coord_equal() + theme_classic()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

#pca <- full_join(pca, zygo_df, by=join_by(ind==Bioinformatic.name))

pve$cumulative <- cumsum(pve$pve)
print(pve$cumulative) #we keep 18 axis

model <- lm(cbind(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18) ~ col, pca)
Manova(model, test.statistic = "Pillai")
model <- lm(cbind(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18) ~ loc, pca)
Manova(model, test.statistic = "Pillai")

# read in data
pca <- read_table2("mortel_GWA.eigenvec", col_names = FALSE)
eigenval <- scan("mortel_GWA.eigenval")

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# sort out the individual location and color
# color
col <- rep(NA, length(pca$ind))
col[grep("blue", pca$ind)] <- "Blue"
col[grep("oran", pca$ind)] <- "Orange"
# location
loc <- rep(NA, length(pca$ind))
loc[grep("CAC", pca$ind)] <- "Cacao"
loc[grep("KAW", pca$ind)] <- "Kaw-Patawa"
loc[grep("PAT", pca$ind)] <- "Kaw-Patawa"
# combine - if you want to plot each in different colours
col_loc <- paste0(col, "_", loc)

# remake data.frame
pca <- as_tibble(data.frame(pca, col, loc, col_loc))

zygo_df <- read.delim('../../name_bioinfo_sample.txt', sep="\t", header= TRUE)
colnames(zygo_df)[2] <- "ind"
zygo_df$Zygosity <- 0
for (i in 1:nrow(pca)){
  if (pca$PC1[i] > 0){
    zygo_df[which(zygo_df$Bioinformatic.name == pca$ind[i]),3] <- "HOM_BL"
  }else if (pca$PC1[i] < (-0.2)){
    zygo_df[which(zygo_df$Bioinformatic.name == pca$ind[i]),3] <- "HOM_OR"
  }else{
    zygo_df[which(zygo_df$Bioinformatic.name == pca$ind[i]),3] <- "HET"
  }
}
write.table(zygo_df, file="../../Zygosity.txt", sep="\t", col.names = TRUE, row.names = FALSE)

# first convert to percentage variance explained
pve <- data.frame(PC = 1:length(eigenval), pve = eigenval/sum(eigenval)*100)
# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_classic()

# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = col, shape = loc)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("blue", "orange"), name="Color") + scale_shape(name="Location")
b <- b + coord_equal() + theme_classic()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

pca <- full_join(pca, zygo_df, by=join_by(ind==Bioinformatic.name))

pve$cumulative <- cumsum(pve$pve)
print(pve$cumulative) #we keep 7 axis
# plot pca with right groups
c <- ggplot(pca, aes(PC1, PC2, col = Zygosity)) + geom_point(size = 3)
c <- c + scale_colour_manual(values = c("goldenrod", "deepskyblue2", "darkorange2"), name="Genotype", labels=c("O/B", "B/B","O/O"))
c <- c + coord_equal() + theme_classic()
c + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

model <- lm(cbind(PC1, PC2, PC3, PC4, PC5, PC6, PC7) ~ Zygosity, pca)
Manova(model, test.statistic = "Pillai")
model <- lm(cbind(PC1, PC2, PC3, PC4, PC5, PC6, PC7) ~ loc, pca)
Manova(model, test.statistic = "Pillai")

## LD decay plotting script

rm(list = ls())

# set path
my_bins <- "mortel_chr1.ld_decay_bins"

# read in data
ld_bins <- read_tsv(my_bins)

# plot LD decay
ggplot(ld_bins, aes(distance/1000, avg_R2)) +
  geom_line(colour="chocolate2") + 
  scale_x_continuous(trans = 'log10', limits=c(NA, 600)) +
  xlab("Distance between SNPs (kbp)") + 
  ylab("Mean Linkage Disequilibrium (LD)") +
  geom_vline(xintercept=c(1e+01,1e+02), linetype="dashed", colour="dodgerblue2", size=1) +
  theme_classic()

