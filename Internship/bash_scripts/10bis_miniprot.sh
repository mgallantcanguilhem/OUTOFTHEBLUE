#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G

##variables
genome_folder="00bis_genomeref"
genome="mortel_blue_ref.softmasked"
genome_or="mortel_orange_ref.softmasked"
file_folder="11_blast_genes"
out_folder="miniprot"
SLURM_DIR="/scratch/mgallantcanguilhem/project_internship"

cd $SLURM_DIR

## Ensure output directory exists
mkdir -p "$file_folder"/"$out_folder"

##load miniprot
source /local/env/envminiprot-0.13.sh

#running different times to find different protein sequence in the genome 

# miniprot -t2 -d "$genome_folder"/"$genome".mpi "$genome_folder"/"$genome".fa
# miniprot -t2 "$genome_folder"/"$genome".mpi "$file_folder"/cortex_prot.fa > "$file_folder"/"$out_folder"/cortex_prot.paf

# miniprot -t2 "$genome_folder"/"$genome".mpi "$file_folder"/OR/OR13_telb_prot.fa > "$file_folder"/OR/OR13_telb_prot.paf
# miniprot -t2 "$genome_folder"/"$genome".mpi "$file_folder"/OR/OR13_telg1_prot.fa > "$file_folder"/OR/OR13_telg1_prot.paf
# miniprot -t2 "$genome_folder"/"$genome".mpi "$file_folder"/OR/OR13_telg2_prot.fa > "$file_folder"/OR/OR13_telg2_prot.paf
# miniprot -t2 "$genome_folder"/"$genome".mpi "$file_folder"/timeless_prot_aligned.fasta > "$file_folder"/"$out_folder"/timeless_blue.paf
# miniprot -t2 "$genome_folder"/"$genome".mpi "$file_folder"/cyr_1_prot_aligned.fasta > "$file_folder"/"$out_folder"/cyr1_blue.paf
# miniprot -t2 "$genome_folder"/"$genome".mpi "$file_folder"/cyr_2_prot_aligned.fasta > "$file_folder"/"$out_folder"/cyr2_blue.paf
# miniprot -t2 "$genome_folder"/"$genome".mpi "$file_folder"/vrille_prot_aligned.fasta > "$file_folder"/"$out_folder"/vrille_blue.paf

# # miniprot -t2 -d "$genome_folder"/"$genome_or".mpi "$genome_folder"/"$genome_or".fa
# miniprot -t2 "$genome_folder"/"$genome_or".mpi "$file_folder"/OR/OR13_telb_prot.fa >> "$file_folder"/OR/OR13_telb_prot.paf
# miniprot -t2 "$genome_folder"/"$genome_or".mpi "$file_folder"/OR/OR13_telg1_prot.fa >> "$file_folder"/OR/OR13_telg1_prot.paf
# miniprot -t2 "$genome_folder"/"$genome_or".mpi "$file_folder"/OR/OR13_telg2_prot.fa >> "$file_folder"/OR/OR13_telg2_prot.paf

# miniprot -t2 -d 16_ragtag_ordered/maniola_jurtina_whole_genome_renamed.mpi 16_ragtag_ordered/maniola_jurtina_whole_genome_renamed.fa
# miniprot -t2 16_ragtag_ordered/maniola_jurtina_whole_genome_renamed.mpi "$file_folder"/optix_prot.fa > "$file_folder"/"$out_folder"/optix_maniola.paf

miniprot -t2 -d 16_ragtag_ordered/spemor_ref.mpi 16_ragtag_ordered/spemor_ref.fna
miniprot -t2 16_ragtag_ordered/spemor_ref.mpi "$file_folder"/optix_prot.fa > "$file_folder"/"$out_folder"/optix_speyeria.paf
