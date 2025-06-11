#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G

##variables
input_BCF="10_VCF_filtered/bgzip/genotyped_mortel_hard_filtered.vcf.bgz"
out_folder="15_het_rate/coverage"
SLURM_DIR="/scratch/mgallantcanguilhem/project_internship/"

cd $SLURM_DIR

##load bcf tools
source /local/env/envbcftools-1.9.sh

# separating my VCF into three VCFs for all genotype
bcftools view $input_BCF -r ptg000030l -S "$out_folder"/HOMBL.txt | sed -E 's/^ptg0*([0-9]+)l/\1/' > "$out_folder"/sampled_HOMBL.vcf
bcftools view $input_BCF -r ptg000030l -S "$out_folder"/HET.txt | sed -E 's/^ptg0*([0-9]+)l/\1/' > "$out_folder"/sampled_HET.vcf
bcftools view $input_BCF -r ptg000030l -S "$out_folder"/HOMOR.txt | sed -E 's/^ptg0*([0-9]+)l/\1/' > "$out_folder"/sampled_HOMOR.vcf

##environment for python
# . $HOME/env/python_winpca/bin/activate

# python "$out_folder"/het_count.py -i "$out_folder"/sampled_HOMBL.vcf -o "$out_folder"/HOMBL_chr30
# python "$out_folder"/het_count.py -i "$out_folder"/sampled_HET.vcf -o "$out_folder"/HET_chr30
# python "$out_folder"/het_count.py -i "$out_folder"/sampled_HOMOR.vcf -o "$out_folder"/HOMOR_chr30

#calculating depth at each sites for each genotype
source /local/env/envvcftools-0.1.16.sh

vcftools --vcf "$out_folder"/sampled_HOMBL.vcf --chr 30 --site-mean-depth --out "$out_folder"/HOMBL_chr30
vcftools --vcf "$out_folder"/sampled_HET.vcf --chr 30 --site-mean-depth --out "$out_folder"/HET_chr30
vcftools --vcf "$out_folder"/sampled_HOMOR.vcf --chr 30 --site-mean-depth --out "$out_folder"/HOMOR_chr30

rm "$out_folder"/sampled_*.vcf #remove the temporary VCF for storage space
