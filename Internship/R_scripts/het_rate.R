setwd('~/Documents/CollègeFr/Stage/Metrics/het_rate')

library(ggplot2)
library(tidyverse)
library(rstatix)
library(ggpubr)

HOMOR <-read.table("HOMOR_chr30_heterozygosity.txt",header=T)
colnames(HOMOR)[2] <- "HOMOR"
HOMBL <-read.table("HOMBL_chr30_heterozygosity.txt",header=T)
colnames(HOMBL)[2] <- "HOMBL"
HET <-read.table("HET_chr30_heterozygosity.txt",header=T)

# # Load PLINK LD file
# ld_forward <- read.table("../11_plink/mortel_chr30.ld", header = TRUE)
# ld_forward <- ld_forward[,c(2,5,7)]
# ld_reversed <- read.table("../11_plink/mortel_chr30_reversed.ld", header = TRUE)
# 
# # Ensure column names match (may vary depending on PLINK version)
# colnames(ld_forward) <- c("POS", "BP_B", "R2")
# colnames(ld_reversed) <- c("POS", "BP_B", "R2")  # Swap BP_A and BP_B
# 
# # Merge both datasets
# ld_combined <- bind_rows(ld_forward, ld_reversed)
# 
# # Compute mean R2 per SNP_A position
# mean_r2 <- ld_combined %>%
#   group_by(POS) %>%
#   summarise(mean_r2 = mean(R2, na.rm = TRUE))
# 
# mean_r2 %>%
#   #filter(POS<1.6 & POS>1.4) %>%
#   mutate(POS=POS/1000000) %>%
#   ggplot(aes(x=POS, y=mean_r2)) +
#   geom_line() +
#   #geom_vline(xintercept = c(1.492798, 1.549730), linewidth=1.3, linetype="longdash") +
#   xlab("Position along Chr 30 (Mb)") +
#   ylab("Linkage Disequilibrium") +
#   theme_classic()

# Load vcftools hap LD file
ld_phased <- read.table("../11_plink/mortel_chr30_10kb.hap.ld", header = TRUE)[, c(2,3,5)]
ld_phased_reversed <- ld_phased[,c(2,1,3)] %>% arrange(POS2)
colnames(ld_phased_reversed)[c(2,1)] <-  colnames(ld_phased)[c(2,1)] 

ld_combined_phased <- bind_rows(ld_phased, ld_phased_reversed) %>% arrange(POS1, POS2)

# Compute mean R2 per POS
mean_r2 <- ld_combined_phased %>%
  group_by(POS1) %>%
  summarise(mean_r2 = mean(R.2, na.rm = TRUE))
colnames(mean_r2)[1] <- "POS"

rate_df <- full_join(HOMOR, HOMBL, by="POS") %>% full_join(HET, by="POS") %>% full_join(mean_r2, by="POS")

rate_df <- rate_df %>% pivot_longer(cols = c(2,3,4), names_to = "Phenotype", values_to = "Heterozygosity_rate") %>% mutate(POS=POS/1000000)

rate_df %>% 
  filter(POS<1.6 & POS>1.4) %>%
  ggplot(aes(x=POS)) + 
  geom_smooth(aes(y=Heterozygosity_rate, col=Phenotype), se=FALSE, span=0.162) +
  geom_smooth(aes(y=mean_r2), col="black", se=FALSE, span=0.162) +
  geom_vline(xintercept = c(1.492798, 1.549730), linewidth=1.3, linetype="longdash") +
  xlab("Position along Chr 30 (Mb)") +
  scale_y_continuous(name ="Rate of heterozygosity", sec.axis = dup_axis(name="Linkage Disequilibrium")) +
  theme_classic()

wilcox <- rate_df %>%
  filter(POS<1.549730 & POS>1.492798) %>%
  wilcox_test(Heterozygosity_rate ~ Phenotype) %>%
  add_significance() %>%
  add_xy_position(x="Phenotype")
wilcox[c(2,3),c(12,13)] <- wilcox[c(3,2),c(12,13)]

rate_df %>% 
  filter(POS<1.549730 & POS>1.492798) %>%
  ggplot(aes(x=Phenotype, y=Heterozygosity_rate)) +
  geom_violin(aes(fill=Phenotype), draw_quantiles = c(0.25,0.5,0.75)) +
  scale_fill_manual(values = c("goldenrod", "deepskyblue2", "darkorange2"), guide="none") +
  geom_jitter(height = 0, width = 0.1) +
  scale_x_discrete(limits=c("HOMBL","HET","HOMOR")) +
  stat_pvalue_manual(wilcox, label= "p.adj", tip.length = 0.01, y.position = c(1.048, 1.048, 1.120), bracket.shorten = 0.05) +
  theme_classic()

#now let's see the coverage across sites
HOMOR_depth <- read.table("HOMOR_chr30.ldepth.mean",header=T)
HOMOR_depth$Geno <- "HOMOR"
HOMBL_depth <-read.table("HOMBL_chr30.ldepth.mean",header=T)
HOMBL_depth$Geno <- "HOMBL"
HET_depth <-read.table("HET_chr30.ldepth.mean",header=T)
HET_depth$Geno <- "HET"

geno_depth_df <- rbind(HOMOR_depth, HOMBL_depth,HET_depth) %>% mutate(POS=POS/1000000)

geno_depth_df %>% 
  filter(POS<1.6 & POS>1.4) %>%
  ggplot(aes(x=POS)) + 
  geom_smooth(aes(y=MEAN_DEPTH, col=Geno), se=FALSE, span=0.2) +
  geom_vline(xintercept = c(1.492798, 1.549730), linewidth=1.3, linetype="longdash") +
  xlab("Position along Chr 30 (Mb)") +
  scale_y_continuous(name ="Coverage") +
  scale_color_manual(values = c("goldenrod", "deepskyblue2", "darkorange2")) +
  theme_classic()

geno_depth_df %>% 
  filter(POS<1.56 & POS>1.48) %>%
  ggplot(aes(x=POS)) + 
  geom_smooth(aes(y=MEAN_DEPTH, col=Geno), se=FALSE, span=0.2) +
  geom_vline(xintercept = c(1.492798, 1.549730), linewidth=1.3, linetype="longdash") +
  xlab("Position along Chr 30 (Mb)") +
  scale_y_continuous(name ="Coverage") +
  scale_color_manual(values = c("goldenrod", "deepskyblue2", "darkorange2")) +
  theme_classic()

#now with the indels
setwd('~/Documents/CollègeFr/Stage/Metrics/het_rate/indel')

HOMOR <-read.table("HOMOR_chr30_heterozygosity.txt",header=T)
colnames(HOMOR)[2] <- "HOMOR"
HOMBL <-read.table("HOMBL_chr30_heterozygosity.txt",header=T)
colnames(HOMBL)[2] <- "HOMBL"
HET <-read.table("HET_chr30_heterozygosity.txt",header=T)

rate_df <- full_join(HOMOR, HOMBL, by="POS") %>% full_join(HET, by="POS") %>% full_join(mean_r2, by="POS")
rate_df <- rate_df %>% pivot_longer(cols = c(2,3,4), names_to = "Phenotype", values_to = "Heterozygosity_rate") %>% mutate(POS=POS/1000000)

rate_df %>% 
  filter(POS<1.6 & POS>1.4) %>%
  ggplot(aes(x=POS)) + 
  geom_smooth(aes(y=Heterozygosity_rate*100, col=Phenotype), se=FALSE, span=0.15) +
  scale_color_manual(values = c("deepskyblue2", "goldenrod", "darkorange2"), name="Genotype", labels=c("B/B", "B/O", "O/O"), limits=c("HOMBL","HET","HOMOR")) +
  geom_smooth(aes(y=mean_r2*100), col="black", se=FALSE, span=0.162) +
  geom_vline(xintercept = c(1.492798, 1.549730), linewidth=1.3, linetype="longdash") +
  xlab("Position along Chr 30 (Mb)") +
  scale_y_continuous(name ="Heterozygosity rate smoothed (%)", sec.axis = sec_axis(~./100 ,name="Linkage Disequilibrium (LD)")) +
  theme_classic()

wilcox <- rate_df %>%
  filter(POS<1.549730 & POS>1.492798) %>%
  wilcox_test(Heterozygosity_rate ~ Phenotype) %>%
  add_significance() %>%
  add_xy_position(x="Phenotype")
wilcox[c(2,3),c(12,13)] <- wilcox[c(3,2),c(12,13)]

rate_df %>% 
  filter(POS<1.549730 & POS>1.492798) %>%
  ggplot(aes(x=Phenotype, y=Heterozygosity_rate*100)) +
  geom_violin(aes(fill=Phenotype), draw_quantiles = c(0.25,0.5,0.75)) +
  scale_fill_manual(values = c("deepskyblue2", "goldenrod", "darkorange2"), name="Genotype", labels=c("B/B", "B/O", "O/O"), limits=c("HOMBL","HET","HOMOR")) +
  geom_jitter(height = 0, width = 0.1) +
  scale_x_discrete(name="Genotype", labels=c("B/B", "B/O", "O/O"), limits=c("HOMBL","HET","HOMOR")) +
  ylab("Heterozygosity rate (%)") +
  stat_pvalue_manual(wilcox, label= "p.adj", tip.length = 0.01, y.position = c(1.048, 1.048, 1.120)*100, bracket.shorten = 0.05) +
  theme_classic()


#now let's see the coverage across sites
HOMOR_depth <- read.table("HOMOR_chr30_indel.ldepth.mean",header=T)
HOMOR_depth$Geno <- "HOMOR"
HOMBL_depth <-read.table("HOMBL_chr30_indel.ldepth.mean",header=T)
HOMBL_depth$Geno <- "HOMBL"
HET_depth <-read.table("HET_chr30_indel.ldepth.mean",header=T)
HET_depth$Geno <- "HET"

geno_depth_df <- rbind(HOMOR_depth, HOMBL_depth,HET_depth) %>% mutate(POS=POS/1000000)

geno_depth_df %>% 
  filter(POS<1.6 & POS>1.4) %>%
  ggplot(aes(x=POS)) + 
  geom_smooth(aes(y=MEAN_DEPTH, col=Geno), se=FALSE, span=0.2) +
  geom_vline(xintercept = c(1.492798, 1.549730), linewidth=1.3, linetype="longdash") +
  xlab("Position along Chr 30 (Mb)") +
  scale_y_continuous(name ="Coverage") +
  scale_color_manual(values = c("goldenrod", "deepskyblue2", "darkorange2")) +
  theme_classic()

geno_depth_df %>% 
  filter(POS<1.56 & POS>1.48) %>%
  ggplot(aes(x=POS)) + 
  geom_smooth(aes(y=MEAN_DEPTH, col=Geno), se=FALSE, span=0.2) +
  geom_vline(xintercept = c(1.492798, 1.549730), linewidth=1.3, linetype="longdash") +
  xlab("Position along Chr 30 (Mb)") +
  scale_y_continuous(name ="Coverage") +
  scale_color_manual(values = c("goldenrod", "deepskyblue2", "darkorange2")) +
  theme_classic()