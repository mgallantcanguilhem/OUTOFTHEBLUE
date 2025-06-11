#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G

##variables
chemin_in1="/projects/outoftheblue/illumina_pitie/PFGS_2024_0948-R1_201224_LOPEZ/fastq_merged_lanes"
chemin_in2="/projects/outoftheblue/illumina_pitie/PFGS_2024_0948-R2_241224_LOPEZ/fastq_merged_lanes"
chemin_out="/projects/outoftheblue/renamed_morpho"

##main code for renaming
tail -n +2 /scratch/$USER/project_internship/correspondance_nom.csv | while IFS=',' read -r col1 col2 #go trough the text of name correspondance
do
    from="$col1" #assign old name to a variable
    to="$col2" #assign new name
    if [[ -e "$chemin_in1/$from" ]]; then #there are two folder with my sequences so I test in one then in another
        cp "$chemin_in1/$from" "$chemin_out/$to" #copy paste the file with a new name in the right place
    else
        cp "$chemin_in2/$from" "$chemin_out/$to"
    fi

done
