#!/bin/bash -l

#SBATCH -A sens2017106
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 3-00:00:00
#SBATCH -J TIDDIT
#SBATCH -C usage_mail

module load bioinfo-tools
module load python
module load pysam

python __main__.py frac_reads /proj/sens2017106/Esmee/bam/2022-07682-02.bam /proj/sens2017106/Esmee/bam/AllNIPT.csv testcount

