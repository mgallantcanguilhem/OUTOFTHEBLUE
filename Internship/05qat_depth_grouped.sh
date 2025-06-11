#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G

# Global variables
input_folder="06_clipped"
output_folder="06ter_depth_grouped"
name_file="names.txt"
SLURM_WD="/scratch/$USER/project_internship/"

cd $SLURM_WD

##environment for samtools
source /local/env/envconda.sh
conda activate $HOME/env/samtools_1.21

#creates the file containing the names of all the clipped aligned reads
> "$input_folder"/"$name_file"

#adds the names
ls -1 "$input_folder"/*.clipped.bam |
while read i
do
    echo "$i" >> "$input_folder"/"$name_file"
done

#now let's run samtools depth on everyone together
samtools depth -f "$input_folder"/"$name_file" | gzip > "$output_folder"/grouped_depth_coverage.gz #get a base depth for all individuals
