#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=assembly
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=32
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4G

hifiasm -o C_sordidulus__KU39667.asm -t 32 C_sordidulus__KU39667.fasta.gz


#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=assembly
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=32
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4G

hifiasm -o C_virens__KU39615.asm -t 32 C_virens__KU39615.fasta.gz
