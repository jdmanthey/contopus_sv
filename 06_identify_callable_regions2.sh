#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=callable2
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G

workdir=/lustre/scratch/jmanthey/18_contopus_sv

input=$( ls ${workdir}/03_callable/*_callable-sorted.bed )

# identify regions of the genome callable from at least 3 (i.e., > 2) samples (max = 1 missing)
bedtools multiinter -i ${input} | \
awk '$4 > 2' | bedtools merge > ${workdir}/03_callable/callable-regions.bed

