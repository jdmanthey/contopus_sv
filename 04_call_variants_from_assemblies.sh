#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=align_call
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=32
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --array=1-4

workdir=/lustre/scratch/jmanthey/18_contopus_sv

basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/haplotypes.txt | tail -n1 )

reference=/home/jmanthey/references/GCA_038380795.1_CAS-CDF_Pyrnan_1.0_genomic.fna

threads=32

# align assemblies to reference genome
minimap2 -cx asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 \
--cs -t ${threads} ${reference} ${workdir}/01_haplotypes/${basename_array}.fa | sort -k6,6 -k8,8n > \
${workdir}/02_aligned/${basename_array}.paf


# call variants from alignments
paftools.js call -L 50000 -s ${basename_array} \
-f ${reference} ${workdir}/02_aligned/${basename_array}.paf | sed 's|1/1|1|g' > \
${workdir}/02_aligned/${basename_array}.vcf


# compress and index vcf 
bgzip -c ${workdir}/02_aligned/${basename_array}.vcf > ${workdir}/02_aligned/${basename_array}.vcf.gz
tabix -p vcf ${workdir}/02_aligned/${basename_array}.vcf.gz








