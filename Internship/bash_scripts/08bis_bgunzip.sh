#!/usr/bin/env bash
#SBATCH -c2
#SBATCH --mem=100G

#variables
gVCF_folder="09_gVCF"
gVCF_file="genotyped_mortel.vcf"
output_folder="09bis_VCF_bgz"
SLURM_WD="/scratch/mgallantcanguilhem/project_internship"

cd $SLURM_WD

##environment
source /local/env/envhtslib-1.6.sh

#bgunzip the file
bgzip -c "$gVCF_folder"/"$gVCF_file" > "$output_folder"/"$gVCF_file".gz

#index the gunzipped VCF file
tabix -p vcf "$output_folder"/"$gVCF_file".gz
