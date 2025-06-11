#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

MEAN_DEPTH="mean_depth.txt"
SLURM_WD="/scratch/mgallantcanguilhem/project_internship/06bis_depth"

cd $SLURM_WD

> $MEAN_DEPTH #creates/cleans the file

# Loop through each .gz file individually
for file in *.gz; do
    # Use zcat to decompress the file and process it
    mdpt=$(zcat "$file" | awk '{ total += $3 } END { if (NR > 0) print total/NR; else print 0 }')
    
    # Write the result to the output file
    echo "${mdpt}\t${file}" >> $MEAN_DEPTH
done
