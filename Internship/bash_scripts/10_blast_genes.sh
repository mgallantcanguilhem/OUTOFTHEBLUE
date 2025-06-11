#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G

##variables
genome_folder="00bis_genomeref"
genome="mortel_blue_ref.hardmasked.fa"
out_folder="11_blast_genes"
SLURM_DIR="/scratch/mgallantcanguilhem/project_internship"

cd $SLURM_DIR

##load blast
source /local/env/envblast-2.12.0.sh

# #create a blast database on the wanted assembly
# makeblastdb -in "$genome_folder"/"$genome" -dbtype nucl -out "$out_folder"/mortel_blue_ref_db

# #run blast
# blastn -query "$out_folder"/cortex.fa -db "$out_folder"/mortel_blue_ref_db -out "$out_folder"/cortex.txt -evalue 1e-5 -outfmt 6
# blastn -query "$out_folder"/optix.fa -db "$out_folder"/mortel_blue_ref_db -out "$out_folder"/optix.txt -evalue 1e-5 -outfmt 6
# blastn -query "$out_folder"/ivory.fa -db "$out_folder"/mortel_blue_ref_db -out "$out_folder"/ivory.txt -evalue 1e-5 -outfmt 6
# blastn -query "$out_folder"/cortex-fly.fa -db "$out_folder"/mortel_blue_ref_db -out "$out_folder"/cortex-fly.txt -evalue 1e-5 -outfmt 6
# blastn -query "$out_folder"/aristaless.fa -db "$out_folder"/mortel_blue_ref_db -out "$out_folder"/aristaless.txt -evalue 1e-5 -outfmt 6
# blastn -query "$out_folder"/aristaless-like.fa -db "$out_folder"/mortel_blue_ref_db -out "$out_folder"/aristaless-like.txt -evalue 1e-5 -outfmt 6
# blastn -query "$out_folder"/yellow.fa -db "$out_folder"/mortel_blue_ref_db -out "$out_folder"/yellow.txt -evalue 1e-5 -outfmt 6
# blastn -query "$out_folder"/tan.fa -db "$out_folder"/mortel_blue_ref_db -out "$out_folder"/tan.txt -evalue 1e-5 -outfmt 6
# blastn -query "$out_folder"/cortex-maniola.fa -db "$out_folder"/mortel_blue_ref_db -out "$out_folder"/cortex-maniola.txt -evalue 1e-5 -outfmt 6
# blastn -query "$out_folder"/cortex-pararge.fa -db "$out_folder"/mortel_blue_ref_db -out "$out_folder"/cortex-pararge.txt -evalue 1e-5 -outfmt 6

# tblastn -query "$out_folder"/wnta_prot.fa -db "$out_folder"/mortel_blue_ref_db -out "$out_folder"/wnta_prot.txt -evalue 1e-5 -outfmt 6
# tblastn -query "$out_folder"/cortex_prot.fa -db "$out_folder"/mortel_blue_ref_db -out "$out_folder"/cortex_prot.txt -evalue 1e-5 -outfmt 6

#makeblastdb -in "$genome_folder"/mortel_orange_ref.softmasked.fa -dbtype nucl -out "$out_folder"/mortel_orange_ref_db
blastn -query "$out_folder"/chr30.fa -db "$out_folder"/mortel_orange_ref_db -out "$out_folder"/chr30.txt -evalue 1e-5 -outfmt 6

# blastn -query "$out_folder"/OR13_telemachus_mrna.fasta -db "$out_folder"/mortel_blue_ref_db -out "$out_folder"/OR13_telemachus_mrna_blue.txt -evalue 1e-5 -outfmt 6
# blastn -query "$out_folder"/OR13_telemachus_mrna.fasta -db "$out_folder"/mortel_orange_ref_db -out "$out_folder"/OR13_telemachus_mrna_orange.txt -evalue 1e-5 -outfmt 6

# blastn -query "$out_folder"/OR/OR13_telb_dna.fasta -db "$out_folder"/mortel_blue_ref_db -out "$out_folder"/OR/OR13_telb.txt -evalue 1e-5 -outfmt 6
# blastn -query "$out_folder"/OR/OR13_telg_1_dna.fasta -db "$out_folder"/mortel_orange_ref_db -out "$out_folder"/OR/OR13_telg1.txt -evalue 1e-5 -outfmt 6
# blastn -query "$out_folder"/OR/OR13_telg_2_dna.fasta -db "$out_folder"/mortel_orange_ref_db -out "$out_folder"/OR/OR13_telg2.txt -evalue 1e-5 -outfmt 6

#makeblastdb -in "$genome_folder"/morrhe_ref.softmasked.fa -dbtype nucl -out "$out_folder"/morrhe_ref_db
#makeblastdb -in "$genome_folder"/morhec_ref.softmasked.fa -dbtype nucl -out "$out_folder"/morhec_ref_db
blastn -query "$out_folder"/chr30.fa -db "$out_folder"/morhec_ref_db -out "$out_folder"/chr30_hec.txt -evalue 1e-5 -outfmt 6
blastn -query "$out_folder"/chr30.fa -db "$out_folder"/morrhe_ref_db -out "$out_folder"/chr30_rhe.txt -evalue 1e-5 -outfmt 6
blastn -query "$out_folder"/optix.fa -db "$out_folder"/morhec_ref_db -out "$out_folder"/optix_hec.txt -evalue 1e-5 -outfmt 6
blastn -query "$out_folder"/optix.fa -db "$out_folder"/morrhe_ref_db -out "$out_folder"/optix_rhe.txt -evalue 1e-5 -outfmt 6
