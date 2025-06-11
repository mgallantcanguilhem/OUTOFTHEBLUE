# Example R Script for simple output plots
# Here, we use pi and dxy output files directly from pixy.
setwd('~/Documents/CollègeFr/Stage/Metrics/10_pixy')

library(ggplot2)
library(tidyverse)

# Provide path to input. Can be pi or Dxy.
# NOTE: this is the only line you should have to edit to run this code:
inp<-read.table("output/pixy_dxy.txt",sep="\t",header=T)
for (i in 1:length(inp$chromosome)){
  inp$chromosome[i] <- substr(inp$chromosome[i], 8, 9)
  if (substr(inp$chromosome[i], 1, 1)=="0"){
    inp$chromosome[i] <- substr(inp$chromosome[i], 2, 2)
  }
}
#write.table(inp, file = "output/long/pixy_fst.txt",sep="\t", row.names=FALSE, quote=FALSE)

# Find the chromosome names and order them: first numerical order, then any non-numerical chromosomes
#   e.g., chr1, chr2, chr22, chrX
for (i in 1:length(inp$chromosome)){
  inp$chromosome[i] <- substr(inp$chromosome[i], 8, 9)
}
chroms <- unique(inp$chromosome)
chrOrder <- sort(chroms)
inp$chrOrder <- factor(inp$chromosome,levels=chrOrder)
inp <- inp %>% filter(chromosome==30)

# Plot pi for each population found in the input file
# Saves a copy of each plot in the working directory
if("avg_pi" %in% colnames(inp)){
  pops <- unique(inp$pop)
  for (p in pops){
    thisPop <- subset(inp, pop == p)
    # Plot stats along all chromosomes:
    popPlot <- ggplot(thisPop, aes(window_pos_1, avg_pi, color=chrOrder)) +
      geom_point()+
      facet_grid(. ~ chrOrder)+
      labs(title=paste("Pi for population", p))+
      labs(x="Position of window start", y="Pi")+
      scale_color_manual(values=rep(c("black","gray"),ceiling((length(chrOrder)/2))))+
      theme_classic()+
      theme(legend.position = "none")
    ggsave(paste("piplot_", p,".png", sep=""), plot = popPlot, device = "png", dpi = 300, width = 14)
  }
} else {
  print("Pi not found in this file")
}

# Plot Dxy for each combination of populations found in the input file
# Saves a copy of each plot in the working directory
if("avg_dxy" %in% colnames(inp)){
  # Get each unique combination of populations
  pops <- unique(inp[c("pop1", "pop2")])
  for (p in 1:nrow(pops)){
    combo <- pops[p,]
    thisPop <- subset(inp, pop1 == combo$pop1[[1]] & pop2 == combo$pop2[[1]])
    # Plot stats along all chromosomes:
    popPlot <- ggplot(thisPop, aes(window_pos_1, avg_dxy, color=chrOrder)) +
      geom_point()+
      facet_grid(. ~ chrOrder)+
      labs(title=paste("Dxy for", combo$pop1[[1]], "&", combo$pop2[[1]]))+
      labs(x="Position of window start", y="Dxy")+
      theme(legend.position = "none")+
      scale_color_manual(values=rep(c("black","gray"),ceiling((length(chrOrder)/2))))+
      theme_classic()+
      theme(legend.position = "none")
    ggsave(paste("dxyplot_", combo$pop1[[1]], "_", combo$pop2[[1]],".png", sep=""), plot = popPlot, device = "png", dpi = 300, width = 14)
  }
} else {
  print("Dxy not found in this file")
}

pixy_to_long <- function(pixy_files){
  
  pixy_df <- list()
  
  for(i in 1:length(pixy_files)){
    
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
    
    if(stat_file_type == "pi"){
      
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value") %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
      
      pixy_df[[i]] <- df
      
      
    } else{
      
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value")
      pixy_df[[i]] <- df
      
    }
    
  }
  
  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
  
}

pixy_folder <- "output/long"
pixy_files <- list.files(pixy_folder, full.names = TRUE)
pixy_df <- pixy_to_long(pixy_files)

# custom labeller for special characters in pi/dxy/fst
pixy_labeller <- as_labeller(c(avg_pi = "pi",
                               avg_dxy = "D[XY]",
                               avg_wc_fst = "F[ST]"),
                             default = label_parsed)

# plotting summary statistics along a single chromosome
chr_num = 30
pixy_df %>%
  filter(chromosome == chr_num) %>%
  filter(statistic %in% c("avg_pi", "avg_dxy", "avg_wc_fst")) %>%
  mutate(chr_position = ((window_pos_1 + window_pos_2)/2)/1000000) %>%
  filter(chr_position<2 & chr_position > 1) %>%
  ggplot(aes(x = chr_position, y = value, color = statistic))+
  geom_line(linewidth = 0.25)+
  facet_grid(statistic ~ .,
             scales = "free_y", switch = "x", space = "free_x",
             labeller = labeller(statistic = pixy_labeller,
                                 value = label_value))+
  xlab(paste0("Position on Chromosome ",chr_num," (Mb)"))+
  ylab("Statistic Value")+
  theme_bw()+
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  scale_color_brewer(palette = "Set1")

# plotting summary statistics across all chromosomes
pixy_df %>%
  mutate(chrom_color_group = case_when(as.numeric(chromosome) %% 2 != 0 ~ "even",
                                       chromosome == "X" ~ "even",
                                       TRUE ~ "odd" )) %>%
  mutate(chromosome = factor(chromosome, levels = c(1:30))) %>%
  filter(statistic %in% c("avg_wc_fst")) %>%
  ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value, color = chrom_color_group))+
  geom_line(linewidth = 0.25)+
  facet_grid(statistic ~ chromosome,
             scales = "free_y", switch = "x", space = "free_x",
             labeller = labeller(statistic = pixy_labeller,
                                 value = label_value))+
  xlab("Chromsome")+
  ylab("Statistic Value")+
  scale_color_manual(values = c("grey50", "black"))+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position ="none")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,NA))

#trying to do a beautiful plot
inp<-read.table("output/pixy_fst.txt",sep="\t",header=T)

# Find the chromosome names and order them: first numerical order, then any non-numerical chromosomes
#   e.g., chr1, chr2, chr22, chrX
for (i in 1:length(inp$chromosome)){
  inp$chromosome[i] <- substr(inp$chromosome[i], 8, 9)
  if (substr(inp$chromosome[i], 1, 1)=="0"){
    inp$chromosome[i] <- substr(inp$chromosome[i], 2, 2)
  }
}

inp$cum_pos <- 0 
inp$cum_pos[1] <- (inp$window_pos_1[1]+inp$window_pos_2[1]+1)/2
for (i in 2:length(inp$cum_pos)){
  if (inp$chromosome[i]==inp$chromosome[i-1]){
    inp$cum_pos[i] <- inp$cum_pos[i-1] + (inp$window_pos_2[i]-inp$window_pos_1[i]+1) + (inp$window_pos_1[i]-inp$window_pos_2[i-1]-1)
  }
  else{
    inp$cum_pos[i] <- inp$cum_pos[i-1] + (inp$window_pos_2[i]-inp$window_pos_1[i]+1) + (inp$window_pos_1[i]-1)
  }
}

chr_label <- c("chr")
chr_x <- c(-3000000)
for (i in unique(inp$chromosome)){
  chr_label <- append(chr_label, i)
  chr_x <- append(chr_x, ((tail(inp[inp$chromosome==i,]$cum_pos,n=1) + inp[inp$chromosome==i,]$cum_pos[1] )/2))
}
chr_df <- data.frame(chr_x, chr_label)

chromosomes <- unique(inp$chromosome)  # Get unique chromosome names in order
chrom_colors <- rep(c("goldenrod", "dodgerblue3"), length.out = length(chromosomes))  # Alternate colors

kb10 <- ggplot(inp, mapping = aes(cum_pos, avg_wc_fst, colour = chromosome)) + 
  geom_line() +
  scale_color_manual(values = setNames(chrom_colors, chromosomes)) +
  scale_x_continuous(name = "Chromosome", breaks = c(1, 2.5, 3.75, 5.5, 7, 8.75, 10.5, 12.5, 14.25, 16, 17.5, 19, 20.5, 22, 23.5, 25.25, 26.75, 28.5, 29.5, 30.25, 31.25, 33, 34.5, 35.5, 37, 39, 39.75, 40.25, 41.25, 42.15, 43.75, 44.75)*10000000, labels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","","","21","22","","24","25","26","","28","","30","32","")) +
  theme_classic() +
  ylab("Fixation index (Fst)") +
  ylim(c(-0.05, 0.45)) +
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position ="none")

#trying to do a beautiful plot
inp<-read.table("output/pixy-window100kb_fst.txt",sep="\t",header=T)

# Find the chromosome names and order them: first numerical order, then any non-numerical chromosomes
#   e.g., chr1, chr2, chr22, chrX
for (i in 1:length(inp$chromosome)){
  inp$chromosome[i] <- substr(inp$chromosome[i], 8, 9)
  if (substr(inp$chromosome[i], 1, 1)=="0"){
    inp$chromosome[i] <- substr(inp$chromosome[i], 2, 2)
  }
}

inp$cum_pos <- 0 
inp$cum_pos[1] <- (inp$window_pos_1[1]+inp$window_pos_2[1]+1)/2
for (i in 2:length(inp$cum_pos)){
  if (inp$chromosome[i]==inp$chromosome[i-1]){
    inp$cum_pos[i] <- inp$cum_pos[i-1] + (inp$window_pos_2[i]-inp$window_pos_1[i]+1) + (inp$window_pos_1[i]-inp$window_pos_2[i-1]-1)
  }
  else{
    inp$cum_pos[i] <- inp$cum_pos[i-1] + (inp$window_pos_2[i]-inp$window_pos_1[i]+1) + (inp$window_pos_1[i]-1)
  }
}

chr_label <- c("chr")
chr_x <- c(-3000000)
for (i in unique(inp$chromosome)){
  chr_label <- append(chr_label, i)
  chr_x <- append(chr_x, ((tail(inp[inp$chromosome==i,]$cum_pos,n=1) + inp[inp$chromosome==i,]$cum_pos[1] )/2))
}
chr_df <- data.frame(chr_x, chr_label)

chromosomes <- unique(inp$chromosome)  # Get unique chromosome names in order
chrom_colors <- rep(c("goldenrod", "dodgerblue3"), length.out = length(chromosomes))  # Alternate colors

kb100 <- ggplot(inp, mapping = aes(cum_pos, avg_wc_fst, colour = chromosome)) + 
  geom_line() +
  scale_color_manual(values = setNames(chrom_colors, chromosomes)) +
  scale_x_continuous(name = "Chromosome", breaks = c(1, 2.5, 3.75, 5.5, 7, 8.75, 10.5, 12.5, 14.25, 16, 17.5, 19, 20.5, 22, 23.5, 25.25, 26.75, 28.5, 29.5, 30.25, 31.25, 33, 34.5, 35.5, 37, 39, 39.75, 40.25, 41.25, 42.15, 43.75, 44.75)*10000000, labels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","","","21","22","","24","25","26","","28","","30","32","")) +
  theme_classic() +
  ylab("Fixation index (Fst)") +
  ylim(c(-0.05, 0.45)) +
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position ="none")

ggarrange(kb10, kb100, align ="hv", ncol = 1, labels = c("A", "B"))

#trying to do a beautiful plot
inp<-read.table("output/pixy-window1kb_fst.txt",sep="\t",header=T)

# Find the chromosome names and order them: first numerical order, then any non-numerical chromosomes
#   e.g., chr1, chr2, chr22, chrX
for (i in 1:length(inp$chromosome)){
  inp$chromosome[i] <- substr(inp$chromosome[i], 8, 9)
  if (substr(inp$chromosome[i], 1, 1)=="0"){
    inp$chromosome[i] <- substr(inp$chromosome[i], 2, 2)
  }
}

inp$cum_pos <- 0 
inp$cum_pos[1] <- (inp$window_pos_1[1]+inp$window_pos_2[1]+1)/2
for (i in 2:length(inp$cum_pos)){
  if (inp$chromosome[i]==inp$chromosome[i-1]){
    inp$cum_pos[i] <- inp$cum_pos[i-1] + (inp$window_pos_2[i]-inp$window_pos_1[i]+1) + (inp$window_pos_1[i]-inp$window_pos_2[i-1]-1)
  }
  else{
    inp$cum_pos[i] <- inp$cum_pos[i-1] + (inp$window_pos_2[i]-inp$window_pos_1[i]+1) + (inp$window_pos_1[i]-1)
  }
}

chr_label <- c("chr")
chr_x <- c(-3000000)
for (i in unique(inp$chromosome)){
  chr_label <- append(chr_label, i)
  chr_x <- append(chr_x, ((tail(inp[inp$chromosome==i,]$cum_pos,n=1) + inp[inp$chromosome==i,]$cum_pos[1] )/2))
}
chr_df <- data.frame(chr_x, chr_label)

chromosomes <- unique(inp$chromosome)  # Get unique chromosome names in order
chrom_colors <- rep(c("black", "gray"), length.out = length(chromosomes))  # Alternate colors

kb1 <- ggplot(inp, mapping = aes(cum_pos, avg_wc_fst, colour = chromosome)) + 
  geom_line() +
  scale_color_manual(values = setNames(chrom_colors, chromosomes)) +
  geom_text(chr_df, mapping = aes(chr_x, rep(0.6, length(chr_x)), label=chr_label), inherit.aes = FALSE, check_overlap = TRUE) +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position ="none")

#combine the FST graphs
ggpubr::ggarrange(kb1, kb10, kb100, nrow=3, labels = c("1kb", "10kb", "100kb"), label.x = c(-0.01, -0.02, -0.02), label.y = 0.2)

#zoom on 30
zoom_30 <- inp %>%
  filter(chromosome=="30" & 1000000<window_pos_1 & window_pos_1<2000000) %>%
  na.omit()

ggplot(zoom_30, mapping = aes(window_pos_1/1000000, avg_wc_fst)) + 
  geom_line(col="brown3", linewidth =1) +
  theme_classic() +
  xlab("Genomic Position on Chromosome 30 (Mb)") +
  ylab("Fixation index (Fst)") +
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position ="none")
