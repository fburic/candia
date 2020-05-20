#!/bin/bash

if [ "$#" -lt 2 ]; then
    echo "[ERROR] Correct syntax:"
	echo "\t" $0 " EXPERIMENT_CONFIG_FILE N_PARALLEL_DECOMP_PER_GPU [SNAKEMAKE_OPTIONS]"
	exit 1
fi

N_CARDS=2


start_partition_on_device () {
    mkdir -p ${TMPDIR}/mps-${device}
    export CUDA_VISIBLE_DEVICES=${device}
    export CUDA_MPS_PIPE_DIRECTORY=${TMPDIR}/mps-${device}/pipe
    export CUDA_MPS_LOG_DIRECTORY=${TMPDIR}/mps-${device}/log

    nvidia-cuda-mps-control -d && printf '\n%s\n' "INFO: MPS daemon started on device ${device}"

    snakemake -p --nolock --forceall --keep-going -j $2 \
      -s parafac_tensorly/decompose_parafac.Snakefile \
      --config npartitions=${N_CARDS} partition=${device} pipe_id=${device} $3 \
      --configfile $1

    echo quit | nvidia-cuda-mps-control && echo "INFO: Stopped MPS daemon on device ${device}"
}


for ((device=0; device<${N_CARDS}; device++)); do

    timestamp=$(date "+%Y%m%d%H%M%S")

    start_partition_on_device $1 $2 $3 \
    > "decompose_partition_${device}_${timestamp}.log" \
    2> "decompose_partition_${device}_${timestamp}.err" &

done
