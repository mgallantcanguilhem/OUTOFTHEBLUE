#!/usr/bin/env bash
#SBATCH -c2
#SBATCH --mem=100G

##environment
source /local/env/envconda.sh
conda activate $HOME/env/bwa-mem2

#variables
genome="mortel_blue_ref.hardmasked.fa"
SLURM_DIR="/scratch/mgallantcanguilhem/project_internship/00bis_genomeref"

##index the genome
bwa-mem2 index "$SLURM_DIR/$genome"
