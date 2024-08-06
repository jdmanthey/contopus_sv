#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=convert
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G

source activate pacbio

bam2fasta -o C_sordidulus__KU39667 m84046_240803_051617_s4.hifi_reads.bc2039.bam 


#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=convert
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G

source activate pacbio

bam2fasta -o C_virens__KU39615 m84046_240803_051617_s4.hifi_reads.bc2040.bam
