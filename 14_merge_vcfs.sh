#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=pg_merge
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=4
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G

source activate bcftools

workdir=/lustre/scratch/jmanthey/18_contopus_sv

input_vcfs=$( ls ${workdir}/06_pangenie_output/*biallelic-filtered2.vcf.gz )

# create a multisample VCF containing all haplotypes
bcftools merge -m id ${input_vcfs} \
> ${workdir}/07_merged_vcf/all_sequence_variants.vcf

