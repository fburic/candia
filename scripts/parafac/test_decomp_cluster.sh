#!/bin/bash

#SBATCH --ntasks 1
#SBATCH --time 00:30:00
#SBATCH --output='decompose-%j.out'
#SBATCH --error='decompose-%j.err'

export CUDA_VISIBLE_DEVICES=0
source parafac_tensorly/start_mps.sh

python parafac_tensorly/decompose_parafac.py

source parafac_tensorly/stop_mps.sh
