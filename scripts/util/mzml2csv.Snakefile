# Note: configfile must contain a list of mzML file paths called 'samples_mzml'
# and a file path called 'swath_windows'

import os.path as path
from glob import glob


input_files = glob(path.join(config['root_dir'], config["samples_mzml"], '*.mzML'))
output_files = [path.join(config['root_dir'], config["samples_csv"],
                          path.split(path.splitext(fname)[0] + '.csv')[1])
                for fname in input_files]


rule all:
    input:
        output_files


rule extract_windows:
    input:
        input_files[0]
    output:
        path.join(config['root_dir'], config["swath_windows"])
    run:
        import convert_mzml2csv as mzml2csv
        mzml2csv.create_isolation_windows_file(str(input), str(output))


rule file_conversion:
    input:
        swath_windows = rules.extract_windows.output,
        mzml_files = path.join(config['root_dir'], config["samples_mzml"], "{sample}.mzML")
    output:
        path.join(config['root_dir'], config["samples_csv"], "{sample}.csv")
    shell:
        "python scripts/util/convert_mzml2csv.py \
        -i {input.mzml_files} \
        -s {input.swath_windows} \
        --min_intensity {config[min_scan_intensity]} \
        -b 5000 \
        -o {output}"
