#!/usr/bin/env bash
#SBATCH --cpus-per-task=6
#SBATCH --mem=400G

# Global variables
GATK="/groups/bipaa/env/gatk4/bin/gatk"
genome_folder="00bis_genomeref"
genome="mortel_blue_ref.hardmasked.fa"
output_folder="09_gVCF/old_gatk"
gVCF_file="genotyped_mortel_old_gatk.vcf"
SLURM_DIR="/scratch/mgallantcanguilhem/project_internship"

cd $SLURM_DIR

# Load needed modules
export JAVA_HOME=/scratch/rpoloni/env/jdk-17.0.6
export PATH=$JAVA_HOME/bin:$PATH
source /local/env/envconda.sh
conda activate /groups/bipaa/env/gatk4


#run GenotypeGVCF and include all non variants sites
$GATK GenotypeGVCFs -R "$genome_folder"/"$genome" -V gendb://08ter_genomicDB_old_gatk --include-non-variant-sites true -O "$output_folder"/"$gVCF_file"

#bgunzip the file
source /local/env/envhtslib-1.6.sh

#compress and index the resulting file
bgzip -c "$output_folder"/"$gVCF_file" > "$output_folder"/"$gVCF_file".gz
tabix -p vcf "$output_folder"/"$gVCF_file".gz
