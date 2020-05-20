import os
from multiprocessing import Pool, cpu_count
from pathlib import Path
import yaml

import pandas as pd


# Set up us some paths
diaumpire_out_dir = Path(config['diaumpire_out_dir'])
msgf_out_dir = diaumpire_out_dir / 'msgf_out'

# Pass the YAML config to the rules that need it by writing it to a temp file
tmp_configfile_copy = f"{os.getenv('TMPDIR')}/config.yaml"
with open(tmp_configfile_copy, 'w') as configfile_copy:
    yaml.dump(config, configfile_copy, default_flow_style = False)


input_files = {
    fname.stem : str(fname) for fname in diaumpire_out_dir.glob('*.mgf')
}
def umpire_sources(wildcards):
    return input_files[wildcards.feature_file]


def load_results(filename):
    return pd.read_table(str(filename)).assign(filename=filename)


rule all:
    input:
        diaumpire_out_dir / 'msgf_ids.csv'


rule msgf:
    """
    Use MS-GF+ to get PSMs from DIA-Umpire pseudo-spectra (feature files)
    """
    input:
        diaumpire_out_dir / '{feature_file}.mgf'
    output:
        msgf_out_dir / '{feature_file}.tsv'
    shell:
        """
        python scripts/util/wrappers.py --process msgf \
        -i {input} -o {output} \
        --config {tmp_configfile_copy}
        """


rule collate_results:
    input:
         expand(rules.msgf.output, feature_file = input_files.keys())
    output:
         diaumpire_out_dir / 'msgf_ids.csv'
    run:
        with Pool(processes=cpu_count()) as pool:
            msgf_results = pool.map(load_results, input)
        msgf_results = pd.concat(msgf_results, ignore_index=True)
        msgf_results.to_csv(str(output), index=False)

