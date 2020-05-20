if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "[ERROR] Don't run $0, source it. Run:  source $0" >&2
    exit 1
fi

pipe_id=$1

mkdir -p ${TMPDIR}/mps-${pipe_id}
export CUDA_MPS_PIPE_DIRECTORY=${TMPDIR}/mps-${pipe_id}/pipe
export CUDA_MPS_LOG_DIRECTORY=${TMPDIR}/mps-${pipe_id}/log
nvidia-cuda-mps-control -d
