#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G

##variables
GATK="/groups/bipaa/env/gatk4/bin/gatk"
genome_folder="00bis_genomeref"
genome="mortel_blue_ref.hardmasked"
seq_folder="06_clipped"
out_folder="07_singlevcf/old_gatk"
SLURM_DIR="/scratch/mgallantcanguilhem/project_internship"

cd $SLURM_DIR

##load module
export JAVA_HOME=/scratch/rpoloni/env/jdk-17.0.6
export PATH=$JAVA_HOME/bin:$PATH
source /local/env/envconda.sh
conda activate /groups/bipaa/env/gatk4

#check genome indexing or index it
if [[ -e ""$genome_folder"/"$genome".fa.fai" ]]; then
    echo "Genome already indexed properly"
else
    #environment for samtools
    source /local/env/envconda.sh
    conda activate $HOME/env/samtools_1.21

    #indexing (again) the genome reference
    samtools faidx "$genome_folder"/"$genome".fa #index the genome

    echo "Genome indexed successfuly"
fi

#check if genome dictionnary exists or creates it
if [[ -e ""$genome_folder"/"$genome".dict" ]]; then
    echo "Genome dictionary already exists"
else
    #indexing (again) the genome reference
    $GATK CreateSequenceDictionary -R "$genome_folder"/"$genome".fa

    echo "Genome dictionary created"
fi

#run haplotype caller in parallel
source /local/env/envparallel-20190122.sh

#creating VCF for each alignment/individual with all sites
ls -1 "$seq_folder"/"$1" |
    parallel -j 10 echo "Creating VFC for {}" \; $GATK HaplotypeCaller -R "$genome_folder"/"$genome".fa -I {} -O "$out_folder"/{/.}.vcf --emit-ref-confidence GVCF --output-mode EMIT_ALL_CONFIDENT_SITES
