#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=vcf_manip
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=32
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4G

source activate bcftools

workdir=/lustre/scratch/jmanthey/18_contopus_sv

input_array=$( ls ${workdir}/03_callable/*txt )

# set alleles outside of callable regions to missing using python script
python3 /home/jmanthey/pangenie_setup/scripts/set-to-missing.py -v ${workdir}/04_pan_vcf/all-haplotypes.vcf \
-m 0.3 -f ${input_array} 2> \
${workdir}/04_pan_vcf/all-haplotypes-callable.log 1> \
${workdir}/04_pan_vcf/all-haplotypes-callable.vcf

# convert haploid VCF into diploid VCF by combining haplotypes of each sample
python3 /home/jmanthey/pangenie_setup/scripts/merge_vcfs.py combine_columns \
-samples ${workdir}/samples.txt -vcf ${workdir}/04_pan_vcf/all-haplotypes-callable.vcf > \
${workdir}/04_pan_vcf/final_callset.vcf


