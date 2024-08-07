#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=align_call
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=20
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=3G

source activate bcftools

workdir=/lustre/scratch/jmanthey/18_contopus_sv

input_vcfs=$( ls ${workdir}/02_aligned/*vcf.gz )

# create a multisample VCF containing all haplotypes
bcftools merge -m none --missing-to-ref ${input_vcfs} | \
python3 /home/jmanthey/pangenie_setup/scripts/assign-variant-ids.py > ${workdir}/04_pan_vcf/all-haplotypes.vcf



