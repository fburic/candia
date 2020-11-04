import os
import numpy as np

# WARNING:
# A wildcard placed at the end of the path string will match
# subdirectories recursively.  Wildcard constraints don't apply here.
#
# So this pattern needs to account for the entire directory hierarchy
# down to known (pattern-free) file names.
slice_dir_pattern = expand("{location}/swath_lower_adjusted={{swath_start}}/rt_window={{time_window}}/slice_tensor.npy",
                           location=config["root_dir"] + '/' + config["slices_location"])

slice_output_pattern = slice_dir_pattern

SWATH_LIST, WINDOW_LIST = glob_wildcards(slice_dir_pattern[0])


def current_partition(file_pattern):
    """
    Split all files matching file_pattern into the config-specified no. of partitions
    and return the partition matching the current config-specified process ID.
    """
    file_list = expand(file_pattern, zip, swath_start=SWATH_LIST, time_window=WINDOW_LIST)
    file_list = np.array(file_list)
    if config.get('subsample', False):
        np.random.seed(42)
        file_list = np.random.choice(file_list, 100, replace=False)

    file_partitions = np.array_split(file_list, config["npartitions"])
    filenames_in_this_partition = [str(filename) for filename in file_partitions[config["partition"]]]

    # Structure from https://stackoverflow.com/a/44590076/10725218
    # Association between output files and source links
    filenames_in_this_partition = {os.path.dirname(fname) + '/parafac.log': fname
                                   for fname in filenames_in_this_partition}
    return filenames_in_this_partition


# Make this association accessible via a function of wildcards
def output_file_to_input_tensor(wildcards):
    partition_filelist = current_partition(slice_dir_pattern)
    return partition_filelist[wildcards.logfile]


rule all:
    input:
        expand('{logfile}', logfile=current_partition(slice_dir_pattern).keys())


rule decompose_slice_gpu:
    # Avoid cyclic rule by making this a parameter
    params:
        tensor = output_file_to_input_tensor,

    output:
        touch('{logfile}')

    shell:
        """
        set +e
        python scripts/parafac/decompose_parafac.py \
        --input_tensor_file {params.tensor} \
        --tensor_properties $(dirname {params.tensor})/tensor_properties.json \
        --init {config[parafac_init]} \
        --min_components {config[parafac_min_comp]} \
        --max_components {config[parafac_max_comp]} \
        --max_iter {config[parafac_max_iter]} \
        --backend {config[parafac_backend]} \
        --avail_ram_gb {config[parafac_avail_ram_gb]} \
        --use_gpu \
        --pipe_id {config[pipe_id]} \
        --skip_completed
        """
