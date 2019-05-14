#!/bin/bash -l

#$ -cwd
#$ -V
#$ -o logs/cluster/
# join stdout and stderr output
#$ -j y
#$ -sync y
#$ -b y
#$ -r y

source /broad/software/scripts/useuse
reuse Bcftools
reuse Samtools
reuse Python-3.4
reuse R-3.5

{exec_job}

