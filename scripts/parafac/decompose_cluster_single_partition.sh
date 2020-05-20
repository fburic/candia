#!/bin/bash

#SBATCH --ntasks 1
#SBATCH --time 54:00:00


if [[ $(hostname) == *"hpc2n.umu.se" ]]; then
  module load icc/2018.3.222-GCC-7.3.0-2.30  impi/2018.3.222
  module load GCC/7.3.0-2.30  CUDA/9.2.88  OpenMPI/3.1.1
  module load TensorFlow/1.10.0-Python-3.6.6

  module load Anaconda3
  source activate parafac3
fi

snakemake --keep-going -s parafac_tensorly/decompose_parafac.Snakefile \
	  --config npartitions=4 partition=3 \
	  --configfile $1
