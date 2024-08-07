#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=callable1
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=20
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-4

workdir=/lustre/scratch/jmanthey/18_contopus_sv

basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/haplotypes.txt | tail -n1 )

reference=/home/jmanthey/references/GCA_038380795.1_CAS-CDF_Pyrnan_1.0_genomic.fna

threads=20

# align assemblies to reference and produce BAM output
minimap2 -ax asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 \
--cs -t ${threads} ${reference} ${workdir}/01_haplotypes/${basename_array}.fa | samtools view -bS | samtools \
sort -o ${workdir}/03_callable/${basename_array}.bam

#index
samtools index ${workdir}/03_callable/${basename_array}.bam

# identify regions of the genome covered by at least one contig from this haplotype
bedtools bamtobed -i ${workdir}/03_callable/${basename_array}.bam | awk '($3-$2) >= 50000' \
| bedtools merge > ${workdir}/03_callable/${basename_array}_covered.bed

# identify regions of the genome that are not covered by more than one contig from this haplotype
bedtools genomecov -bga -ibam ${workdir}/03_callable/${basename_array}.bam | \
awk '$4 < 2' | bedtools merge > ${workdir}/03_callable/${basename_array}_notduplicated.bed

# identify regions of the genome with exactly one contig from this haplotype mapping
# using bedtools interect and the previous two bed files
bedtools intersect -a ${workdir}/03_callable/${basename_array}_covered.bed \
-b ${workdir}/03_callable/${basename_array}_notduplicated.bed > \
${workdir}/03_callable/${basename_array}_callable.bed

# sort the bed file
bedtools sort -i ${workdir}/03_callable/${basename_array}_callable.bed > \
${workdir}/03_callable/${basename_array}_callable-sorted.bed




