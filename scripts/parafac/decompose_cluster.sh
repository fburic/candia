#!/bin/bash

#SBATCH --ntasks 1
#SBATCH --time 168:00:00
#SBATCH --output='decompose-%j.out'
#SBATCH --error='decompose-%j.err'


if [[ $(hostname) == *"hpc2n.umu.se" ]]; then
  module load icc/2018.3.222-GCC-7.3.0-2.30  impi/2018.3.222
  module load GCC/7.3.0-2.30  CUDA/9.2.88  OpenMPI/3.1.1
  module load TensorFlow/1.10.0-Python-3.6.6

  module load Anaconda3
  source activate parafac3
fi

if [[ $(hostname) == "vera"* ]]; then
    N_PAR_DECOMP=8
fi

CUDA_VISIBLE_DEVICES=0

source parafac_tensorly/start_mps.sh

# Meant to be run as an array job.
# The number of partitions is that of the number of array jobs.
# Each job in the array is assigned a data partition using the array task ID.
snakemake -j ${N_PAR_DECOMP} --keep-going -s parafac_tensorly/decompose_parafac.Snakefile \
  --config npartitions=${SLURM_ARRAY_TASK_COUNT} partition=${SLURM_ARRAY_TASK_ID} device=0 \
  --configfile $1

source parafac_tensorly/stop_mps.sh
