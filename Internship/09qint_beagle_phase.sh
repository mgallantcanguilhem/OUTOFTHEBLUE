#!/bin/bash
#SBATCH -c2
#SBATCH --mem=50G

#set variables
BEAGLE_DIR="/scratch/rpoloni/env"
INPUT_BCF="bgzip/genotyped_mortel_hard_filtered.vcf.bgz"
out_folder="phased_vcf"
SLURM_WD="/scratch/mgallantcanguilhem/project_internship/10_VCF_filtered"

cd $SLURM_WD
mkdir -p $out_folder

# # Load environments
# export JAVA_HOME=/scratch/rpoloni/env/jdk-17.0.6
# export PATH=$JAVA_HOME/bin:$PATH

# source /local/env/envbcftools-1.9.sh

# # #run beagle genome-wide
# # java -jar "$BEAGLE_DIR"/beagle.22Jul22.46e.jar gt="$INPUT_VCF" out="$OUTDIR"/05_phased impute=true nthreads=2 window=1 overlap=0.1 gp=false

# #run on chromosome 30
# bcftools view -s ^mortel_GUY_PAT_male_oran_10_L1 -r ptg000030l "$INPUT_BCF" | sed -E 's/^ptg0*([0-9]+)l/\1/' > "$out_folder"/chr30.vcf
# java -jar "$BEAGLE_DIR"/beagle.22Jul22.46e.jar gt="$out_folder"/chr30.vcf out="$out_folder"/chr30_phased impute=true nthreads=4 window=1 overlap=0.1 gp=false
# rm "$out_folder"/chr30.vcf

source /local/env/envvcftools-0.1.16.sh
#calculate the linkage disequilibrium with a 10kb window on phased data
vcftools --gzvcf "$out_folder"/chr30_phased.vcf.gz --hap-r2 --ld-window-bp 10000 --out "$out_folder"/mortel_chr30_10kb
