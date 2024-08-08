#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=pg_biallelic
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=4
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-35

source activate bcftools

workdir=/lustre/scratch/jmanthey/18_contopus_sv

basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames.txt | tail -n1 )

# convert pangenie output to biallelic to unmerge any merged sequence variants
cat ${workdir}/06_pangenie_output/${basename_array}_genotyping.vcf | python3 \
/home/jmanthey/pangenie/scripts/convert-to-biallelic.py ${workdir}/04_pan_vcf/final_callset.vcf > \
${workdir}/06_pangenie_output/${basename_array}_biallelic.vcf

# gzip and index
bgzip ${workdir}/06_pangenie_output/${basename_array}_biallelic.vcf

tabix -p vcf ${workdir}/06_pangenie_output/${basename_array}_biallelic.vcf.gz

# filter based on genotype quality per individual (min. GQ = 30)
bcftools view -i 'FMT/GQ>29' \
${workdir}/06_pangenie_output/${basename_array}_biallelic.vcf.gz > \
${workdir}/06_pangenie_output/${basename_array}_biallelic-filtered.vcf

# gzip and index
bgzip ${workdir}/06_pangenie_output/${basename_array}_biallelic-filtered.vcf

tabix -p vcf ${workdir}/06_pangenie_output/${basename_array}_biallelic-filtered.vcf.gz

# make a tmp file for reheader
echo ${basename_array} > ${basename_array}.txt

# reheader the vcf (individual name is just sample, so need to rename)
bcftools reheader -s ${basename_array}.txt -o \
${workdir}/06_pangenie_output/${basename_array}_biallelic-filtered2.vcf.gz \
${workdir}/06_pangenie_output/${basename_array}_biallelic-filtered.vcf.gz

# index
tabix -p vcf ${workdir}/06_pangenie_output/${basename_array}_biallelic-filtered2.vcf.gz

# remove unneeded file
rm ${basename_array}.txt


