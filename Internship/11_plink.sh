#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G

##variables
input_VCF="/scratch/mgallantcanguilhem/project_internship/10_VCF_filtered/genotyped_mortel_hard_filtered.vcf.gz"
input_BCF="/scratch/mgallantcanguilhem/project_internship/10_VCF_filtered/bgzip/genotyped_mortel_hard_filtered.vcf.bgz"
folder_bfile="input_file"
result_folder="GWAS_result"
PCA_folder="PCA_result"
LD_folder="LD_result"
SLURM_DIR="/scratch/mgallantcanguilhem/project_internship/12_plink_GWAS"

cd $SLURM_DIR

##load bcf tools and p link
source /local/env/envconda.sh
conda activate $HOME/env/plink_1.9

source /local/env/envbcftools-1.9.sh

#change contig name from ptg000030l to 30 for instance
zcat $input_VCF | sed -E 's/^ptg0*([0-9]+)l/\1/' > "$folder_bfile"/converted.vcf
bcftools view $input_BCF -r ptg000030l | sed -E 's/^ptg0*([0-9]+)l/\1/' > "$folder_bfile"/sampled_converted.vcf

#calulating LD with another method than the one previously used 1kb window
source /local/env/envvcftools-0.1.16.sh
vcftools --vcf "$folder_bfile"/sampled_converted.vcf --geno-r2 --ld-window-bp 1000 --out "$LD_folder"/mortel_chr30

#creating binary files for plink
plink --vcf "$folder_bfile"/sampled_converted.vcf --make-bed --double-id --set-missing-var-ids @:# --chr-set 68 no-x no-y no-xy no-mt --out "$folder_bfile"/mortel_GWA
plink --vcf "$folder_bfile"/converted.vcf --make-bed --double-id --set-missing-var-ids @:# --chr-set 68 no-x no-y no-xy no-mt --out "$folder_bfile"/mortel
rm "$folder_bfile"/*converted.vcf #remove for storage space the vcf with chromosome name converted

#removing missing sex info to all males
awk '{$5=1; print}' "$folder_bfile"/mortel.fam > temp.fam && mv temp.fam "$folder_bfile"/mortel.fam

#running p link with different types of GWAS, quali and quanti data...
plink --bfile "$folder_bfile"/mortel --chr-set 67 no-x no-y no-xy no-mt --pheno phenotype_qual.txt --assoc --out "$result_folder"/gwas_results

plink --bfile "$folder_bfile"/mortel --chr-set 67 no-x no-y no-xy no-mt --pheno phenotype_brightness.txt --linear --out "$result_folder"/gwas_results_brightness
#plink --bfile "$folder_bfile"/mortel --chr-set 67 no-x no-y no-xy no-mt --keep ind_blue.txt --pheno phenotype_hue_blue.txt --assoc --pfilter 1e-2 --out "$result_folder"/gwas_results_hue

plink --bfile "$folder_bfile"/mortel --chr-set 67 no-x no-y no-xy no-mt --keep ind_blue.txt --pheno phenotype_blue_PC1.txt --linear --out "$result_folder"/gwas_results_blue_PC1

# plink --bfile "$folder_bfile"/mortel --chr-set 67 no-x no-y no-xy no-mt --pheno phenotype_PC1.txt --linear --pfilter 1e-2 --out "$result_folder"/gwas_results_quanti

#plink --bfile "$folder_bfile"/mortel --chr-set 67 no-x no-y no-xy no-mt --pheno phenotype_PC1.txt --linear --out "$result_folder"/gwas_results_quanti_all_snps

# plink --bfile "$folder_bfile"/mortel --chr-set 67 no-x no-y no-xy no-mt --indep-pairwise 50 10 0.1 --out "$folder_bfile"/mortel
# plink --bfile "$folder_bfile"/mortel --chr-set 67 no-x no-y no-xy no-mt --extract "$folder_bfile"/mortel.prune.in --pca --out "$PCA_folder"/mortel

#plink --bfile "$folder_bfile"/mortel_GWA --chr-set 67 no-x no-y no-xy no-mt --indep-pairwise 50 10 0.1 --out "$folder_bfile"/mortel_GWA
#plink --bfile "$folder_bfile"/mortel_GWA --chr-set 67 no-x no-y no-xy no-mt --pca --out "$PCA_folder"/mortel_GWA

# #just quick calulation of LD decay
#plink --bfile "$folder_bfile"/mortel --chr-set 67 no-x no-y no-xy no-mt --chr 1 --thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --out "$LD_folder"/mortel_chr1
#source /local/env/envpython-2.7.15.sh
#python ld_decay_calc.py -i "$LD_folder"/mortel_chr1.ld.gz -o "$LD_folder"/mortel_chr1

#back to gwas
#plink --bfile "$folder_bfile"/mortel --chr-set 67 no-x no-y no-xy no-mt --pheno phenotype_PC1.txt --linear --covar "$PCA_folder"/mortel.eigenvec --pfilter 1e-2 --out "$result_folder"/gwas_results_quanti_pca

#plink --bfile "$folder_bfile"/mortel --chr-set 67 no-x no-y no-xy no-mt --pheno phenotype_PC2.txt --linear --pfilter 1e-2 --out "$result_folder"/gwas_results_quanti2

#calulating the number of independent SNPs to adjust the bonferroni correction
# zcat $input_VCF | sed -E 's/^ptg0*([0-9]+)l/\1/' > "$folder_bfile"/converted.vcf
# plink --vcf "$folder_bfile"/converted.vcf --make-bed --double-id --set-missing-var-ids @:# --chr-set 68 no-x no-y no-xy no-mt --fill-missing-a2 --out "$folder_bfile"/mortel_filled
# rm "$folder_bfile"/converted.vcf
# awk '{$5=1; print}' "$folder_bfile"/mortel_filled.fam > temp.fam && mv temp.fam "$folder_bfile"/mortel_filled.fam
# plink --bfile "$folder_bfile"/mortel_filled --chr-set 67 no-x no-y no-xy no-mt --maf --recode A --out meff_matrix
#source /local/env/envpython-2.7.15.sh
#python trans_meff.py

#back to GWAS
#plink --bfile "$folder_bfile"/mortel --chr-set 67 no-x no-y no-xy no-mt --pheno phenotype_PC1_new.txt --linear --pfilter 1e-2 --out "$result_folder"/gwas_results_quanti_new
#plink --bfile "$folder_bfile"/mortel --chr-set 67 no-x no-y no-xy no-mt --pheno phenotype_PC2_new.txt --linear --pfilter 1e-2 --out "$result_folder"/gwas_results_quanti2_new

# plink --bfile "$folder_bfile"/mortel --chr-set 67 no-x no-y no-xy no-mt --pheno phenotype_PC1_hind.txt --linear --pfilter 1e-2 --out "$result_folder"/gwas_results_quanti_hind
# plink --bfile "$folder_bfile"/mortel --chr-set 67 no-x no-y no-xy no-mt --pheno phenotype_PC2_hind.txt --linear --pfilter 1e-2 --out "$result_folder"/gwas_results_quanti2_hind

plink --bfile "$folder_bfile"/mortel --chr-set 67 no-x no-y no-xy no-mt --pheno phenotype_PC1_all.txt --linear --covar "$PCA_folder"/mortel.eigenvec --out "$result_folder"/gwas_results_quanti_all
plink --bfile "$folder_bfile"/mortel --chr-set 67 no-x no-y no-xy no-mt --pheno phenotype_PC2_all.txt --linear --out "$result_folder"/gwas_results_quanti2_all

#another LD method
# plink --bfile "$folder_bfile"/mortel --chr-set 67 no-x no-y no-xy no-mt --chr 30 -r2 --ld-window 20 --ld-window-kb 10 --ld-window-r2 0 --out "$LD_folder"/mortel_chr30
# awk '{print $2, $5, $7}' "$LD_folder"/mortel_chr30.ld | awk '{print $2, $1, $3}' > "$LD_folder"/mortel_chr30_reversed.ld
