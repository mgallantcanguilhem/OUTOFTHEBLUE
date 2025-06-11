#!/usr/bin/env bash
#SBATCH -c2
#SBATCH --mem=100G

#variables
gVCF_folder="09_gVCF/old_gatk"
gVCF_file="genotyped_mortel_old_gatk.vcf"
metric_folder="09ter_filter_metrics/old_gatk"
subset_VCF="mortel_subset_old.vcf"
OUT="mortel_subset_old"
SLURM_WD="/scratch/mgallantcanguilhem/project_internship"

cd $SLURM_WD

source /local/env/envconda.sh
conda activate $HOME/env/vcflib_1.0.3

#randomly subsample the vcf
vcfrandomsample "$gVCF_folder"/"$gVCF_file" -r 0.012 > "$metric_folder"/"$subset_VCF"


cd $metric_folder

# compress vcf
source /local/env/envhtslib-1.6.sh
bgzip $subset_VCF

source /local/env/envbcftools-1.9.sh
# index vcf
bcftools index "$subset_VCF".gz

#calculate our metrics
source /local/env/envvcftools-0.1.16.sh

vcftools --gzvcf "$subset_VCF".gz --freq2 --out $OUT --max-alleles 2 #MAF
vcftools --gzvcf "$subset_VCF".gz --depth --out $OUT #coverage/ind
vcftools --gzvcf "$subset_VCF".gz --site-mean-depth --out $OUT #coverage per site
vcftools --gzvcf "$subset_VCF".gz --site-quality --out $OUT #quality
vcftools --gzvcf "$subset_VCF".gz --missing-indv --out $OUT #% missing
vcftools --gzvcf "$subset_VCF".gz --missing-site --out $OUT #% missing/site
vcftools --gzvcf "$subset_VCF".gz --het --out $OUT
