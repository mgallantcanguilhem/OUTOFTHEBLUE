#!/usr/bin/env bash
#SBATCH -c2
#SBATCH --mem=100G

#variables
output_folder="bgzip"
SLURM_WD="/scratch/mgallantcanguilhem/project_internship/10_VCF_filtered"

cd $SLURM_WD

##environment
source /local/env/envhtslib-1.6.sh

## Run bgzip and tabix on all filtered files
ls -1 genotyped_mortel_no_indel_filtered.vcf.gz |
while read i
do
    file=$(basename "$i" .vcf.gz)

    echo "Zipping $file"

    zcat "$file".vcf.gz | iconv -f ISO-8859-1 -t UTF-8 | bgzip -c > "$output_folder"/"$file".vcf.bgz
    tabix -p vcf "$output_folder"/"$file".vcf.bgz

    echo "Done for $file"
done