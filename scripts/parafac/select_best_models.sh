#!/bin/bash

#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --time 00:30:00

set -eo pipefail

time Rscript scripts/parafac/select_best_models.R $*
