#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=sniffles2
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=32
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4G

workdir=/lustre/scratch/jmanthey/18_contopus_sv

threads=32

# activate environment
source activate sniffles

# initial multi-sample genotype calling
sniffles --input \
${workdir}/08_sniffles/C_sordidulus__KU39667.snf \
${workdir}/08_sniffles/C_virens__KU39615.snf \
--vcf ${workdir}/08_sniffles/sniffles_genotypes.vcf

