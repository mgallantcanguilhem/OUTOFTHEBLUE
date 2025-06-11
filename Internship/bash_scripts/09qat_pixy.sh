#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G

# Global variables
VCFbzip_folder="10_VCF_filtered/concatenated"
VCFbzip_file="genotyped_mortel_concat.vcf.bgz"
output_folder="10bis_pixy"
SLURM_WD="/scratch/mgallantcanguilhem/project_internship"

cd $SLURM_WD

#load pixy environment
source /local/env/envconda.sh
conda activate $HOME/env/pixy_1.2.11

#running pixy to calculate pi, fst and dxy
pixy --stats pi fst dxy --vcf "$VCFbzip_folder"/"$VCFbzip_file" --populations population.txt --window_size 100000 --n_cores 2 --output_folder $output_folder --outpout_prefix pixy-window100kb
