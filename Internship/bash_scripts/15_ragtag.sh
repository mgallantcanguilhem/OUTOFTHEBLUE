#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G

##variables
genome_folder="/scratch/mgallantcanguilhem/project_internship/00bis_genomeref"
genome_ref="mortel_blue_ref.softmasked.fa"
SLURM_DIR="/scratch/mgallantcanguilhem/project_internship/16_ragtag_ordered"

cd $SLURM_DIR

##load minimap and samtools
source /local/env/envconda.sh
conda activate $HOME/env/ragtag_2.1.0

#reorder ptg according to chromosome assembly of another species

#ragtag.py scaffold $genome_folder/$genome_ref helmel_ref.fna -o Heliconius -t 2 -w
#ragtag.py scaffold $genome_folder/$genome_ref spemor_ref.fna -o Speyeria -t 2 -w

ragtag.py scaffold maniola_jurtina_whole_genome_renamed.fa $genome_folder/$genome_ref -o CHR_ordering -t 2 -w
