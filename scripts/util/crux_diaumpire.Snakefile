import inspect
import sys
from pathlib import Path

_this_filename = inspect.getframeinfo(inspect.currentframe()).filename
_this_path = Path(_this_filename).parent.resolve()
sys.path.append(str(_this_path.parent))

from util import wrappers


INPUT_MAPS, = glob_wildcards(str(Path(config['samples_mzml']) / '{swath_map}.mzML'))


rule all:
    input:
        str(Path(config['diaumpire_out_dir']) / 'diaumpire_crux_psms.csv')

        
rule diaumpire:
    input:
        expand(str(Path(config['samples_mzml']) / '{swath_map}.mzML'), swath_map=INPUT_MAPS)
    output:
        expand(str(Path(config['diaumpire_out_dir']) / '{swath_map}_Q1.mgf'), swath_map=INPUT_MAPS),
        expand(str(Path(config['diaumpire_out_dir']) / '{swath_map}_Q2.mgf'), swath_map=INPUT_MAPS),
        expand(str(Path(config['diaumpire_out_dir']) / '{swath_map}_Q3.mgf'), swath_map=INPUT_MAPS)
    run:
        # TODO: parallelize
        for input_file in input:
            wrappers.run_diaumpire(input_file, config['diaumpire_out_dir'], config['diaumpire_paramfile'])


rule crux:
    input:
        rules.diaumpire.output
    output:
        expand(str(Path(config['diaumpire_out_dir']) / '{swath_map}_Q1_crux.tsv'), swath_map=INPUT_MAPS),
        expand(str(Path(config['diaumpire_out_dir']) / '{swath_map}_Q2_crux.tsv'), swath_map=INPUT_MAPS),
        expand(str(Path(config['diaumpire_out_dir']) / '{swath_map}_Q3_crux.tsv'), swath_map=INPUT_MAPS)
    run:
        # TODO: parallelize
        for input_file in input:
            wrappers.run_crux(input_file, output_psm_filename=None, config=config)


rule collate_results:
    input:
        rules.crux.output
    output:
        str(Path(config['diaumpire_out_dir']) / 'diaumpire_crux_psms.csv')
    run:
        import pandas as pd
        import re
        psm_tables = []
        for percolator_file in input:
            quality_num = re.findall('_Q\d+', percolator_file)[0].split('_')[-1][-1]
            psms = pd.read_csv(percolator_file, sep='\t')
            psms = psms.assign(file = percolator_file)
            psms = psms.assign(quality = quality_num)
            psm_tables.append(psms)
        psm_tables = pd.concat(psm_tables, ignore_index=True)
        psm_tables.to_csv(str(output), index=False)

