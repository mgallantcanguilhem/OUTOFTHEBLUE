#!/usr/bin/env bash
#SBATCH --cpus-per-task=5
#SBATCH --mem=200G

# Global variables
GATK="/groups/bipaa/env/gatk4/bin/gatk"
genome_folder="00bis_genomeref"
genome="mortel_blue_ref.hardmasked.fa"
VCF_folder="07_singlevcf/old_gatk"
output_folder="/scratch/mgallantcanguilhem/project_internship/08ter_genomicDB_old_gatk"
info_folder="/scratch/mgallantcanguilhem/project_internship/08bis_infoGDB"
SLURM_DIR="/scratch/mgallantcanguilhem/project_internship"

cd $SLURM_DIR

# Load needed modules
export JAVA_HOME=/scratch/rpoloni/env/jdk-17.0.6
export PATH=$JAVA_HOME/bin:$PATH
source /local/env/envconda.sh
conda activate /groups/bipaa/env/gatk4

# ###creates a list of vcf files (output from gatk_hapcaller) and file identifiers

# ls -1 "$VCF_folder"/*vcf > "$info_folder"/sample_list
# orig_list="$info_folder"/sample_list
# modified_list="$info_folder"/"modified_sample_list"
# search_string="07_singlevcf/"
# replace_string=""
# suffix=".clipped.vcf"

# #loops in the file and removes the search_string and vcf suffix
# while IFS= read -r line; do
#     if [[ $line == $search_string* ]]; then
#         line=${line//$search_string/$replace_string}
#         line=${line%$suffix}
#     fi
#     echo "$line" >> "$modified_list"
# done < "$orig_list"

# echo "File names modified successfully."

# #merges the two lists in a tab-delimited file
# paste -d$'\t' "$modified_list" "$orig_list" > "$info_folder"/vcf_list

# echo "vcf list created"
# rm $orig_list
# rm $modified_list

# ## creates a list of scaffolds for genomicsdbimport (yes, you have to specify this even for a genome-wide run)
# grep ptg "$genome_folder"/"$genome" > $info_folder/intervals.list
# sed -i 's/>//g' $info_folder/intervals.list

###creates a genomic database with all vcfs for genotyping
echo "making GenomicsDB"

$GATK GenomicsDBImport \
--sample-name-map $info_folder/vcf_list \
--genomicsdb-workspace-path $output_folder \
-L $info_folder/intervals.list
