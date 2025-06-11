#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G

##variables
genome_folder="00bis_genomeref"
genome="mortel_blue_ref.hardmasked.fa"
seq_folder="02_cleanedseq"
out_folder="03_mapping"
SLURM_DIR="/scratch/mgallantcanguilhem/project_internship"

cd $SLURM_DIR

#loop for cleaning data progressively
for read1 in $(find $seq_folder -name "$1" | sort); 
do
    #find paired read each time
    read2=$(echo "$read1" | perl -pe 's/_R1/_R2/') #find corresponding read2
    echo "Aligning file $read1 $read2"

    #getting just file names and not the path
    name1=$(basename "$read1")
    name2=$(basename "$read2")

    #getting the read group for the reads
    ID="$(echo $name1 | cut -d'_' -f 6-7)" #find ID of reads + lane of Illumina
    RG="@RG\tID:"$ID"_L002\tSM:$ID\tPL:Illumina" #defining the read group

    ##environment for bwa
    source /local/env/envconda.sh
    conda activate $HOME/env/bwa-mem2

    ##map on the indexed genome
    bwa-mem2 mem -R "$RG" "$genome_folder/$genome" "$seq_folder/$name1" "$seq_folder/$name2" > "$out_folder"/"${name1%_cleaned_R1.fastq.gz}.bam"

    ##environment for samtools
    source /local/env/envconda.sh
    conda activate $HOME/env/samtools_1.21

    #cd "$SLURM_DIR/$out_folder"

    #filter alignement over q10
    samtools view -Sb -q 10 "$out_folder"/"${name1%_cleaned_R1.fastq.gz}.bam" | samtools sort -o "$out_folder"/"${name1%_cleaned_R1.fastq.gz}.sorted.bam"

    #index the output alignment
    samtools index "$out_folder"/"${name1%_cleaned_R1.fastq.gz}.sorted.bam"

    #remove unsorted bam (save storage space)
    rm "$out_folder"/"${name1%_cleaned_R1.fastq.gz}.bam"
done
