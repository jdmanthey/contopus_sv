#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=pg_genotype
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=64
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-35

workdir=/lustre/scratch/jmanthey/18_contopus_sv

basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames.txt | tail -n1 )

threads=64

# run pangenie for each sample
PanGenie -f index \
-i <(zcat ${workdir}/05_fastq/${basename_array}_R1.fastq.gz ${workdir}/05_fastq/${basename_array}_R2.fastq.gz) \
-o ${workdir}/06_pangenie_output/${basename_array} -j ${threads} -t ${threads}

