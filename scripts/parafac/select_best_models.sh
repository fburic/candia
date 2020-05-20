#!/bin/bash

#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --time 00:30:00

set -eo pipefail

if [[ $(hostname) == "vera"* ]]; then
    source deactivate
    module load iccifort/2018.3.222-GCC-7.3.0-2.30 impi/2018.3.222 R/3.5.1
fi

time Rscript parafac_tensorly/select_best_models.R $*
