#!/bin/bash

#SBATCH --ntasks 32
#SBATCH --time 24:00:00
#SBATCH --output='diann-%j.out'
#SBATCH --error='diann-%j.err'

snakemake -p --forceall -s quantification/diann.Snakefile $*

