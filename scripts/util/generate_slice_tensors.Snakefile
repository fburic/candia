# Note: configfile must contain a list of mzML file paths called 'slices_location'

import os
import socket

# WARNING:
# A wildcard placed at the end of the path string will match
# subdirectories recursively.  Wildcard constraints don't apply here
# so this can't be constrained.
#
# So this pattern needs to account for the entire directory hierarchy
# down to known (pattern-free) file names.
slice_part_pattern = expand("{location}/swath_lower_adjusted={{swath_start}}/rt_window={{time_window}}/part-{{part}}.csv",
                           location=config["root_dir"] + '/' + config["slices_location"])

SWATH_LIST, WINDOW_LIST, _ = glob_wildcards(slice_part_pattern[0])

slice_dir_pattern = os.path.join(*slice_part_pattern[0].split("/")[0:-1])

# Dummy filename (not a parameter for the script):
# It makes it easier for Snakemake to link rule output with input.
# The name matches one output from the script to avoid junk files.

output_file_name = "slice_tensor.npy"
slice_output_pattern = os.path.join(slice_dir_pattern, output_file_name)


localrules: all

rule all:
    input:
        expand(slice_output_pattern,
               zip, swath_start=SWATH_LIST, time_window=WINDOW_LIST)

rule generate_slice_tensor:
    input:
        slice_dir_pattern

    output:
        touch(slice_output_pattern)

    shell:
        "python scripts/util/generate_slice_tensor.py --slice_dir {input} --mz_tolerance {config[mass_tol_ppm]}"

