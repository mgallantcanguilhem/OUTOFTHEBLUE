#!/usr/bin/env bash
#SBATCH -c2
#SBATCH --mem=100G

#variables
WINPCA="$HOME/env/winpca/winpca"
input_VCF="/scratch/mgallantcanguilhem/project_internship/10_VCF_filtered/genotyped_mortel_hard_filtered.vcf.gz"
biall_VCF="/scratch/mgallantcanguilhem/project_internship/10_VCF_filtered/genotyped_mortel_hard_filtered_biallelic.vcf"
SLURM_WD="/scratch/mgallantcanguilhem/project_internship/13_winpca"

cd $SLURM_WD

##environment
#source /local/env/envbcftools-1.9.sh

#bcftools view -m2 -M2 -v snps $input_VCF -Oz | gunzip -c | iconv -f ISO-8859-1 -t UTF-8 > "$biall_VCF" #putting the vcf in the right format with only biallelic sites

##environment
. $HOME/env/python_winpca/bin/activate

# Run Windowed PCA on all chromosome 30
$WINPCA pca "$biall_VCF" ptg000030l:1-16143720 ptg30_all -w 100000 -i 10000

# Make a plot of principal component 1 and color by population (here morph)
$WINPCA chromplot ptg30_all ptg000030l:1-16143720 -m population.txt -g color
