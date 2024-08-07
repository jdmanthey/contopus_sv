#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=site_ID
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --array=1-4

workdir=/lustre/scratch/jmanthey/18_contopus_sv

basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/haplotypes.txt | tail -n1 )

# extract the 3rd column (variant ID) from the VCF for all callable sites per individual
bedtools intersect -a ${workdir}/04_pan_vcf/all-haplotypes.vcf \
-b ${workdir}/03_callable/${basename_array}_callable.bed \
-wa -f 1.0 | cut -f3 > ${workdir}/03_callable/${basename_array}.txt
