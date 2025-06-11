#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G

##variables
seq_folder="05_dedup"
out_folder="06_clipped"
SLURM_DIR="/scratch/mgallantcanguilhem/project_internship"

cd $SLURM_DIR

for file in $(find "$seq_folder" -name "*.bam" | sort); #going trough my files
do
    echo treating "$file"

    ##environment for bamutil
    source /local/env/envconda.sh
    conda activate $HOME/env/bamutil_1.0.15

    # Clip overlap
    bam clipOverlap \
        --in "$file" \
        --out "$out_folder"/$(basename "$file" .sorted.dedup.bam).temp.bam \
        --unmapped --storeOrig OC --stats --poolSize 2000000

    ##environment for samtools
    source /local/env/envconda.sh
    conda activate $HOME/env/samtools_1.21

    # Remove the reads that became unmapped in the clipping
    samtools view -hb -F 4 \
        "$out_folder"/$(basename "$file" .sorted.dedup.bam).temp.bam  \
        > "$out_folder"/$(basename "$file" .sorted.dedup.bam).clipped.bam

    # Index bam
    samtools index "$out_folder"/$(basename "$file" .sorted.dedup.bam).clipped.bam

    # Cleanup
    rm "$out_folder"/$(basename "$file" .sorted.dedup.bam).temp.bam

done
