#!/usr/bin/env bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=300G

# Global variables
ALIGNEDFOLDER="03_mapping"
METRICSFOLDER="05bis_metrics_dedup"
DEDUPFOLDER="05_dedup"
SLURM_SUBMIT_DIR="/scratch/$USER/project_internship/"

cd $SLURM_SUBMIT_DIR

##environment for parallel
source /local/env/envparallel-20190122.sh

# Remove duplicates from bam alignments in parallel
mkdir tmp 2>/dev/null #making a folder to store temporary files

#removing duplicates in my alignments
ls -1 "$ALIGNEDFOLDER"/*.sorted.bam | tail -n +2 |
    parallel -j 5 echo "Deduplicating sample {}" \; source /local/env/envconda.sh \; conda activate $HOME/env/picard_3.3.0 \; picard MarkDuplicates INPUT={} OUTPUT="$DEDUPFOLDER"/{/.}.dedup.bam METRICS_FILE="$METRICSFOLDER"/{/.}.metrics.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true MAX_FILE_HANDLES=100 TMP_DIR="./tmp" 
rm -r ./tmp #removing the folder
