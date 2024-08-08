#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G

source activate bcftools

# define main working directory
workdir=/lustre/scratch/jmanthey/18_contopus_sv/07_merged_vcf

# filter based on missing data (max missing = 3/35) and MAC > 0 
vcftools --vcf ${workdir}/all_sequence_variants.vcf --max-missing 0.9 \
--max-maf 0.49 --mac 1 --recode --recode-INFO-all \
--out ${workdir}/all_sequence_variants_filtered

# filter based on missing data (max missing = 3/35) and MAC > 9 
vcftools --vcf ${workdir}/all_sequence_variants.vcf --max-missing 0.9 \
--max-maf 0.49 --mac 9 --recode --recode-INFO-all \
--out ${workdir}/all_sequence_variants_filtered_mac10

# filter based on missing data (no missing) and MAC > 9 
vcftools --vcf ${workdir}/all_sequence_variants.vcf --max-missing 1.0 \
--max-maf 0.49 --mac 9 --recode --recode-INFO-all \
--out ${workdir}/all_sequence_variants_filtered_mac10_complete
