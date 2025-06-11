#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G

# Global variables
script_folder="scripts"
data_folder="06ter_depth_grouped"
SLURM_WD="/scratch/$USER/project_internship/"

cd $SLURM_WD

source /local/env/envr-4.3.3.sh

gzip -dc "$data_folder"/grouped_depth_coverage.gz > "$data_folder"/grouped_depth_coverage.txt #making information temporary accesible for my R script

Rscript "$script_folder"/coverage_sliding_wind_cluster.R "$data_folder"/grouped_depth_coverage.txt coverage_sliding_output #running my R script for computing depth

rm "$data_folder"/grouped_depth_coverage.txt #remove temporary file
