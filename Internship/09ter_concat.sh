#!/usr/bin/env bash
#SBATCH -c2
#SBATCH --mem=100G

#variables
file1="bgzip/genotyped_mortel_hard_filtered.vcf.bgz"
file2="bgzip/genotyped_mortel_invariants.vcf.bgz"
output_folder="concatenated"
SLURM_WD="/scratch/mgallantcanguilhem/project_internship/10_VCF_filtered"

cd $SLURM_WD

##environment
source /local/env/envbcftools-1.9.sh
source /local/env/envhtslib-1.6.sh

# combine the two VCFs (variant and invariant sites) using bcftools concat
bcftools concat --allow-overlaps $file1 $file2 -O v | bgzip -c > "$output_folder"/genotyped_mortel_concat.vcf.bgz
tabix -p vcf "$output_folder"/genotyped_mortel_concat.vcf.bgz #indexing the resulting compressed vcf
