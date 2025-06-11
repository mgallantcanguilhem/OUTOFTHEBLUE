#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G

# Global variables
VCFzip_folder="09bis_VCF_bgz"
VCFzip_file="genotyped_mortel.vcf.gz"
output_folder="10_VCF_filtered/"
output_file1="genotyped_mortel_hard_filtered.vcf.gz"
output_file2="genotyped_mortel_no_indel_filtered.vcf.gz"
output_file3="genotyped_mortel_no_MAF_filtered.vcf.gz"
output_file4="genotyped_mortel_no_MAF_lowmiss_filtered.vcf.gz"
output_file5="genotyped_mortel_no_MAF_noqual_filtered.vcf.gz"
output_file6="genotyped_mortel_invariants.vcf.gz"
SLURM_DIR="/scratch/mgallantcanguilhem/project_internship"

cd $SLURM_DIR

#loading needed package
source /local/env/envvcftools-0.1.16.sh

# set filters
MAF=0.1
#careful for MISS you wanna encode 1-(minimum missing freq wanted) !
MISS=0.85
QUAL=25
MIN_DEPTH=7
MAX_DEPTH=40

#> "$output_folder"/num_position.txt #create a file to log the number of sites left after filtering

#echo "$VCFzip_file" >> "$output_folder"/num_position.txt
#gzip -dc "$VCFzip_folder"/"$VCFzip_file" | grep -v "#" | wc -l >> "$output_folder"/num_position.txt

# perform the filtering with vcftools (different type)

#hard filtering
#vcftools --gzvcf "$VCFzip_folder"/"$VCFzip_file" --remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > "$output_folder"/"$output_file1"

#echo "$output_file1" >> "$output_folder"/num_position.txt
#gzip -dc "$output_folder"/"$output_file1" | grep -v "#" | wc -l >> "$output_folder"/num_position.txt

#keep indels
#vcftools --gzvcf "$VCFzip_folder"/"$VCFzip_file" --maf $MAF --max-missing $MISS --minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > "$output_folder"/"$output_file2"

#echo "$output_file2" >> "$output_folder"/num_position.txt
#gzip -dc "$output_folder"/"$output_file2" | grep -v "#" | wc -l >> "$output_folder"/num_position.txt

#keep all non variant site (no MAF filtering)
#vcftools --gzvcf "$VCFzip_folder"/"$VCFzip_file" --remove-indels --max-missing $MISS --minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > "$output_folder"/"$output_file3"

#echo "$output_file3" >> "$output_folder"/num_position.txt
#gzip -dc "$output_folder"/"$output_file3" | grep -v "#" | wc -l >> "$output_folder"/num_position.txt

#inversing missing
#vcftools --gzvcf "$VCFzip_folder"/"$VCFzip_file" --remove-indels --max-missing 0.15 --minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > "$output_folder"/"$output_file4"

#echo "$output_file4" >> "$output_folder"/num_position.txt
#gzip -dc "$output_folder"/"$output_file4" | grep -v "#" | wc -l >> "$output_folder"/num_position.txt

#inversing missing
#vcftools --gzvcf "$VCFzip_folder"/"$VCFzip_file" --remove-indels --max-missing 0.15 --max-meanDP $MAX_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > "$output_folder"/"$output_file5"

#echo "$output_file5" >> "$output_folder"/num_position.txt
#gzip -dc "$output_folder"/"$output_file5" | grep -v "#" | wc -l >> "$output_folder"/num_position.txt

#only invariant sites
vcftools --gzvcf "$VCFzip_folder"/"$VCFzip_file" --max-maf 0 --recode --stdout | gzip -c > "$output_folder"/"$output_file6"

echo "$output_file6" >> "$output_folder"/num_position.txt
gzip -dc "$output_folder"/"$output_file6" | grep -v "#" | wc -l >> "$output_folder"/num_position.txt
