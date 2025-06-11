#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G

# Global variables
GENOMEFOLDER="00bis_genomeref"
GENOME="mortel_blue_ref.hardmasked.fa"
ALIGNEDFOLDER="03_mapping"
METRICSFOLDER="04_metrics_alignment"
SLURM_SUBMIT_DIR="/scratch/$USER/project_internship/"

cd $SLURM_SUBMIT_DIR

##environment for picard
source /local/env/envconda.sh
conda activate $HOME/env/picard_3.3.0

# Run Picard Tools to get some metrics on data & alignments
ls -1 "$ALIGNEDFOLDER"/*.sorted.bam | tail -n +2 | #going trough my files
while read i
do
    file=$(basename "$i") #getting file name

    echo "Computing alignment metrics for $file"

    #running picard to collect metrics
    picard CollectAlignmentSummaryMetrics \
        R="$GENOMEFOLDER"/"$GENOME" \
        I="$ALIGNEDFOLDER"/"$file" \
        O="$METRICSFOLDER"/"$(basename "$file" .sorted.bam)"_alignment_metrics.txt

    echo "Computing insert size metrics for $file"

    picard CollectInsertSizeMetrics \
        I="$ALIGNEDFOLDER"/"$file" \
        OUTPUT="$METRICSFOLDER"/"$(basename "$file" .sorted.bam)"_insert_size_metrics.txt \
        HISTOGRAM_FILE="$METRICSFOLDER"/"$(basename "$file" .sorted.bam)"_insert_size_histogram.pdf

    echo "Done for $file"
done
