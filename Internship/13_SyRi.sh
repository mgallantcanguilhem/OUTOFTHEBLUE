#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G

##variables
ptg_folder="/scratch/mgallantcanguilhem/project_internship/00ter_ptg30"
inv_seq_folder="/scratch/mgallantcanguilhem/project_internship/00qat_inv_seq"
SLURM_DIR="/scratch/mgallantcanguilhem/project_internship/14_SyRI/individual_test"

cd $SLURM_DIR

#load minimap and samtools
source /local/env/envconda.sh
conda activate $HOME/env/samtools_1.21

source /local/env/envminimap2-2.15.sh

#Create BAM file with all the alignments I want on a graph

#blue vs orange
minimap2 -ax asm5 --eqx "$ptg_folder"/ptg30_blue_zoom.fa "$ptg_folder"/ptg27_orange_zoom.fa | samtools sort -O BAM -o "$ptg_folder"/ptg30_zoom.bam
samtools index "$ptg_folder"/ptg30_zoom.bam

#blue vs niepelti vs orange
minimap2 -ax asm5 --eqx "$ptg_folder"/ptg30_blue_zoom_hardmasked.fa "$ptg_folder"/ptg06_nie_zoom_hardmasked.fa | samtools sort -O BAM -o "$ptg_folder"/ptg30_nie-blue_zoom.bam
samtools index "$ptg_folder"/ptg30_nie-blue_zoom.bam
minimap2 -ax asm5 --eqx "$ptg_folder"/ptg06_nie_zoom_hardmasked.fa "$ptg_folder"/ptg27_orange_zoom.fa | samtools sort -O BAM -o "$ptg_folder"/ptg30_nie-orange_zoom.bam
samtools index "$ptg_folder"/ptg30_nie-orange_zoom.bam

minimap2 -ax asm5 --eqx "$ptg_folder"/ptg30_blue_zoom_hardmasked.fa "$ptg_folder"/ptg25_the_zoom_hardmasked.fa | samtools sort -O BAM -o "$ptg_folder"/ptg30_the-blue_zoom.bam
samtools index "$ptg_folder"/ptg30_the-blue_zoom.bam
minimap2 -ax asm5 --eqx "$ptg_folder"/ptg25_the_zoom_hardmasked.fa "$ptg_folder"/ptg27_orange_zoom.fa | samtools sort -O BAM -o "$ptg_folder"/ptg30_the-orange_zoom.bam
samtools index "$ptg_folder"/ptg30_the-orange_zoom.bam

minimap2 -ax asm5 --eqx "$ptg_folder"/ptg30_blue_zoom_hardmasked.fa "$ptg_folder"/ptg17_hec_zoom_hardmasked.fa | samtools sort -O BAM -o "$ptg_folder"/ptg30_hec-blue_zoom.bam
samtools index "$ptg_folder"/ptg30_hec-blue_zoom.bam
minimap2 -ax asm5 --eqx "$ptg_folder"/ptg17_hec_zoom_hardmasked.fa "$ptg_folder"/ptg27_orange_zoom.fa | samtools sort -O BAM -o "$ptg_folder"/ptg30_hec-orange_zoom.bam
samtools index "$ptg_folder"/ptg30_hec-orange_zoom.bam

minimap2 -ax asm5 --eqx "$ptg_folder"/ptg30_blue_zoom_hardmasked.fa "$ptg_folder"/ptg43_rhe_zoom_hardmasked.fa | samtools sort -O BAM -o "$ptg_folder"/ptg30_rhe-blue_zoom.bam
samtools index "$ptg_folder"/ptg30_rhe-blue_zoom.bam
minimap2 -ax asm5 --eqx "$ptg_folder"/ptg43_rhe_zoom_hardmasked.fa "$ptg_folder"/ptg27_orange_zoom.fa | samtools sort -O BAM -o "$ptg_folder"/ptg30_rhe-orange_zoom.bam
samtools index "$ptg_folder"/ptg30_rhe-orange_zoom.bam

minimap2 -ax asm5 --eqx "$ptg_folder"/ptg30_blue_zoom_hardmasked.fa "$ptg_folder"/ptg26_cyp_zoom_hardmasked.fa | samtools sort -O BAM -o "$ptg_folder"/ptg30_cyp-blue_zoom.bam
samtools index "$ptg_folder"/ptg30_cyp-blue_zoom.bam
minimap2 -ax asm5 --eqx "$ptg_folder"/ptg26_cyp_zoom_hardmasked.fa "$ptg_folder"/ptg27_orange_zoom.fa | samtools sort -O BAM -o "$ptg_folder"/ptg30_cyp-orange_zoom.bam
samtools index "$ptg_folder"/ptg30_cyp-orange_zoom.bam

#load SyRI
conda activate $HOME/env/syri_1.6

#run syri with the alignments that I want to find synteny or structural variants
#blue vs orange
syri -c "$ptg_folder"/ptg30_zoom.bam -r "$ptg_folder"/ptg30_blue_zoom.fa -q "$ptg_folder"/ptg27_orange_zoom.fa -k -F B --prefix ptg30_zoom

#blue vs niepelti vs orange
syri -c "$ptg_folder"/ptg30_nie-blue_zoom.bam -r "$ptg_folder"/ptg30_blue_zoom_hardmasked.fa -q "$ptg_folder"/ptg06_nie_zoom_hardmasked.fa -k -F B --prefix ptg30_nie-blue_zoom
syri -c "$ptg_folder"/ptg30_nie-orange_zoom.bam -r "$ptg_folder"/ptg06_nie_zoom_hardmasked.fa -q "$ptg_folder"/ptg27_orange_zoom.fa -k -F B --prefix ptg30_nie-orange_zoom

syri -c "$ptg_folder"/ptg30_the-blue_zoom.bam -r "$ptg_folder"/ptg30_blue_zoom_hardmasked.fa -q "$ptg_folder"/ptg25_the_zoom_hardmasked.fa -k -F B --prefix ptg30_the-blue_zoom
syri -c "$ptg_folder"/ptg30_the-orange_zoom.bam -r "$ptg_folder"/ptg25_the_zoom_hardmasked.fa -q "$ptg_folder"/ptg27_orange_zoom.fa -k -F B --prefix ptg30_the-orange_zoom

syri -c "$ptg_folder"/ptg30_hec-blue_zoom.bam -r "$ptg_folder"/ptg30_blue_zoom_hardmasked.fa -q "$ptg_folder"/ptg17_hec_zoom_hardmasked.fa -k -F B --prefix ptg30_hec-blue_zoom
syri -c "$ptg_folder"/ptg30_hec-orange_zoom.bam -r "$ptg_folder"/ptg17_hec_zoom_hardmasked.fa -q "$ptg_folder"/ptg27_orange_zoom.fa -k -F B --prefix ptg30_hec-orange_zoom

syri -c "$ptg_folder"/ptg30_rhe-blue_zoom.bam -r "$ptg_folder"/ptg30_blue_zoom_hardmasked.fa -q "$ptg_folder"/ptg43_rhe_zoom_hardmasked.fa -k -F B --prefix ptg30_rhe-blue_zoom
syri -c "$ptg_folder"/ptg30_rhe-orange_zoom.bam -r "$ptg_folder"/ptg43_rhe_zoom_hardmasked.fa -q "$ptg_folder"/ptg27_orange_zoom.fa -k -F B --prefix ptg30_rhe-orange_zoom

syri -c "$ptg_folder"/ptg30_cyp-blue_zoom.bam -r "$ptg_folder"/ptg30_blue_zoom_hardmasked.fa -q "$ptg_folder"/ptg26_cyp_zoom_hardmasked.fa -k -F B --prefix ptg30_cyp-blue_zoom
syri -c "$ptg_folder"/ptg30_cyp-orange_zoom.bam -r "$ptg_folder"/ptg26_cyp_zoom_hardmasked.fa -q "$ptg_folder"/ptg27_orange_zoom.fa -k -F B --prefix ptg30_cyp-orange_zoom


#load Plotsr
source /local/env/envplotsr-0.5.3.sh

#plotting the results with color information on the genomes file
plotsr --sr ptg30_zoomsyri.out --genomes ptg30_zoom.txt -W 5 -H 3 -o ptg30_zoom_plot.png

plotsr --sr ptg30_nie-blue_zoomsyri.out --sr ptg30_nie-orange_zoomsyri.out --genomes ptg30_nie.txt -W 6 -H 3 -o ptg30_nie_plot.png
plotsr --sr ptg30_the-blue_zoomsyri.out --sr ptg30_the-orange_zoomsyri.out --genomes ptg30_the.txt -W 6 -H 3 -o ptg30_the_plot.png
plotsr --sr ptg30_hec-blue_zoomsyri.out --sr ptg30_hec-orange_zoomsyri.out --genomes ptg30_hec.txt -W 6 -H 3 -o ptg30_hec_plot.png
plotsr --sr ptg30_rhe-blue_zoomsyri.out --sr ptg30_rhe-orange_zoomsyri.out --genomes ptg30_rhe.txt -W 6 -H 3 -o ptg30_rhe_plot.png
plotsr --sr ptg30_cyp-blue_zoomsyri.out --sr ptg30_cyp-orange_zoomsyri.out --genomes ptg30_cyp.txt -W 6 -H 3 -o ptg30_cyp_plot.png
