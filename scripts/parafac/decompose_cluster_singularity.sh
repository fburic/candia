#!/bin/bash

#SBATCH --ntasks 8
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

# A single GPU is used, on which parallel decompositions will be run

if [[ $(hostname) == "vera"* ]]; then
    N_PAR_DECOMP=6
fi

# Working on a shared machine, so careful to use a unique MPS pipe
pipe_id=$(uuidgen)

singularity exec --nv parafac_tensorly/parafac.simg /bin/bash -c \
  "CUDA_VISIBLE_DEVICES=0 \
  && source parafac_tensorly/start_mps.sh ${pipe_id} \
  && snakemake -j ${N_PAR_DECOMP} --nolock --forceall --keep-going \
  -s parafac_tensorly/decompose_parafac.Snakefile \
  --config npartitions=${SLURM_ARRAY_TASK_COUNT} partition=${SLURM_ARRAY_TASK_ID} pipe_id=${pipe_id} \
  --configfile $1 \
  && source parafac_tensorly/stop_mps.sh"
