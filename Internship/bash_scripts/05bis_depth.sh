#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G

# Global variables
input_folder="06_clipped"
output_folder="06bis_depth"
SLURM_WD="/scratch/$USER/project_internship/"

cd $SLURM_WD

##environment for picard
source /local/env/envconda.sh
conda activate $HOME/env/samtools_1.21

## Run Samtools to get coverage
ls -1 "$input_folder"/*.clipped.bam | #run trough my files
while read i
do
    file=$(basename "$i" .clipped.bam) #get the file name

    echo "Computing coverage for $file" 

    samtools depth "$input_folder"/"$file".clipped.bam | gzip > "$output_folder"/"$file".coverage.gz #retrieving the depth for each file

    echo "Done for $file"
done
