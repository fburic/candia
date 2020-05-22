#!/bin/bash

#SBATCH --ntasks 1
#SBATCH --time 168:00:00
#SBATCH --output='decompose-%j.out'
#SBATCH --error='decompose-%j.err'

if [ "$#" -lt 2 ]; then
    echo "[ERROR] Correct syntax:"
	echo $0 " EXPERIMENT_CONFIG_FILE N_PARALLEL_DECOMP_PER_GPU"
	exit 1
fi

# A single GPU is used for each array job, on which parallel decompositions will be run
N_PAR_DECOMP=$2

CUDA_VISIBLE_DEVICES=0

source scripts/parafac/start_mps.sh

# Meant to be run as an array job.
# The number of partitions is that of the number of array jobs.
# Each job in the array is assigned a data partition using the array task ID.
snakemake -j ${N_PAR_DECOMP} --keep-going -s scripts/parafac/decompose_parafac.Snakefile \
  --config npartitions=${SLURM_ARRAY_TASK_COUNT} partition=${SLURM_ARRAY_TASK_ID} device=0 \
  --configfile $1

source scripts/parafac/stop_mps.sh
