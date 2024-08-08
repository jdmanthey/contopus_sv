#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=pg_index
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=32
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=3G

workdir=/lustre/scratch/jmanthey/18_contopus_sv

reference=/home/jmanthey/references/GCA_038380795.1_CAS-CDF_Pyrnan_1.0_genomic.fna

threads=32

# pangenie preprocessing to make index
PanGenie-index -v ${workdir}/04_pan_vcf/final-graph.vcf -r ${reference} -t ${threads} -o index

