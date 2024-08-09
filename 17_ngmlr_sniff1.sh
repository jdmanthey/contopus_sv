#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=align_call
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=64
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-2

workdir=/lustre/scratch/jmanthey/18_contopus_sv

reference=/home/jmanthey/references/GCA_038380795.1_CAS-CDF_Pyrnan_1.0_genomic.fna

pacbio_base=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames_pacbio.txt | tail -n1 )

threads=64

# align pacbio reads to reference
ngmlr -t ${threads} --bam-fix -r ${reference} \
-q ${workdir}/00_data/${pacbio_base}.fasta.gz  \
-o ${workdir}/08_sniffles/${pacbio_base}.sam

# activate environment
source activate bcftools

# convert sam to bam
samtools view -b -S -o ${workdir}/08_sniffles/${pacbio_base}.bam \
${workdir}/08_sniffles/${pacbio_base}.sam

# remove sam file
rm ${workdir}/08_sniffles/${pacbio_base}.sam

# clean up the bam file
java -jar /home/jmanthey/picard.jar CleanSam \
I=${workdir}/08_sniffles/${pacbio_base}.bam \
O=${workdir}/08_sniffles/${pacbio_base}_cleaned.bam

# remove the raw bam
rm ${workdir}/08_sniffles/${pacbio_base}.bam

# sort the cleaned bam file
java -jar /home/jmanthey/picard.jar SortSam \
I=${workdir}/08_sniffles/${pacbio_base}_cleaned.bam \
O=${workdir}/08_sniffles/${pacbio_base}_cleaned_sorted.bam SORT_ORDER=coordinate

# remove the cleaned bam file
rm ${workdir}/08_sniffles/${pacbio_base}_cleaned.bam

# index the final bam file
samtools index ${workdir}/08_sniffles/${pacbio_base}_cleaned_sorted.bam

# activate environment
source activate sniffles

# initial genotyping of SVs using sniffles
sniffles --input ${workdir}/08_sniffles/${pacbio_base}_cleaned_sorted.bam \
--snf ${workdir}/08_sniffles/${pacbio_base}.snf --threads ${threads} --non-germline






