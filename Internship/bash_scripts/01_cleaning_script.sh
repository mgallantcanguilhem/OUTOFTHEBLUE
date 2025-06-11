#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G

##environment
source /local/env/envfastp-0.20.0.sh

##variables
chemin_in="/projects/outoftheblue/renamed_morpho"
chemin_out="/scratch/$USER/project_internship/02_cleanedseq"

#get input files for paired read R1
reads1=$(find $chemin_in -name '*R1.fastq.gz' | sort)

#loop for cleaning data progressively
for read1 in $reads1; do #going trough the reads one by one
    
    #find paired read each time
    ID="${read1: -17:-15}" #find ID of read1
    read2=$(find $chemin_in -name "*$ID*R2.fastq.gz") #find corresponding read2

    name_out1="$chemin_out/${read1: -42:-12}_cleaned_R1" #changes the name of read 1
    name_out2="$chemin_out/${read1: -42:-12}_cleaned_R2" #changes the name of read 2
    name_report="fastp_$ID" #stock the name of the read

    #run fastp to clean the reads
    fastp -i "$read1" -I "$read2" -o "$name_out1.fastq.gz" -O "$name_out2.fastq.gz" -h "$name_report.html" -j "$name_report.json"
done
