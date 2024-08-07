#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=graph_construct
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=32
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4G

source activate bcftools

workdir=/lustre/scratch/jmanthey/18_contopus_sv

reference=/home/jmanthey/references/GCA_038380795.1_CAS-CDF_Pyrnan_1.0_genomic.fna

# construct graph using merge_vcfs python script
python3 /home/jmanthey/pangenie_setup/scripts/merge_vcfs.py merge \
-vcf ${workdir}/04_pan_vcf/final_callset.vcf \
-r ${reference} -ploidy 2 2> ${workdir}/04_pan_vcf/graph-filtered-tmp.log 1> \
${workdir}/04_pan_vcf/final-graph-tmp.vcf

# gzip and index
bgzip ${workdir}/04_pan_vcf/final-graph-tmp.vcf
tabix -p vcf ${workdir}/04_pan_vcf/final-graph-tmp.vcf.gz

# normalize the indels in the VCF
bcftools norm -d all -f ${reference} ${workdir}/04_pan_vcf/final-graph-tmp.vcf.gz > \
${workdir}/04_pan_vcf/final-graph.vcf

