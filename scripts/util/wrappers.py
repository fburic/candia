import argparse
import inspect
import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
import sys
import yaml

import coloredlogs

_this_filename = inspect.getframeinfo(inspect.currentframe()).filename
_this_path = Path(_this_filename).parent.resolve()
sys.path.append(str(_this_path.parent))

from ..util import msproc

logging.basicConfig(format=msproc.LOG_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)
coloredlogs.install(fmt=msproc.LOG_FORMAT, level='INFO', logger=logger)


def run_diaumpire(input_filename, output_dir, param_filename):
    """
    Assumes the DIAUMPIRE_PATH env var is set with the location of DIA-Umpire
    """
    DIAUMPIRE_PATH = os.getenv('DIAUMPIRE_PATH', '')
    cmd = (f'java -jar -Xmx8G {DIAUMPIRE_PATH}/DIA_Umpire_SE.jar '
           f'{input_filename} {param_filename}')
    _run_cmd(cmd, 'Running DIA-Umpire:')

    # Wildcards are not used to allow for parallel calls of this function
    # (each call should be concerned with ony its files)
    input_path = Path(input_filename)
    output_filename_pattern = str(input_path.parent / input_path.stem) + '_Q{}.mgf'
    for qual_num in range(1, 4):
        try:
            shutil.move(output_filename_pattern.format(qual_num), output_dir)
        except FileNotFoundError as e:
            logger.warning(f'File {output_filename_pattern.format(qual_num)} was not found. Skipping.')
    logger.info(f'Output files moved to {output_dir}')


def run_crux(input_filename, output_psm_filename, config):
    """
    Get Comet hits against target and decoy databases, then use Percolator
    to separate the high-confidence PSMs.
    The target and decoy hits are created as temporary files only..

    :param input_filename: Scan file for which to get PSMs
    :param output_psm_filename: High-confidence PSMs file.
           By setting it to None, the output filename will be at the same location as the input
           and have the same name, except with a '_crux.tsv' suffix instead of '.mgf'
    :param config: Dictionary-like structure holding necessary params for the Crux tools
    """
    output_targets_filename = tempfile.mkstemp()[1]
    output_decoys_filename = tempfile.mkstemp()[1]
    run_comet(input_filename, output_targets_filename, config['database'],
              mass_tol_ppm=config['mass_tol_ppm'])
    run_comet(input_filename, output_decoys_filename, config['decoy_database'],
              mass_tol_ppm=config['mass_tol_ppm'])

    if output_psm_filename is None:
        output_psm_filename = str(Path(input_filename).parent / Path(input_filename).stem) + '_crux.tsv'
    run_percolator(output_targets_filename, output_decoys_filename, output_psm_filename,
                   config['percolator_fdr'], config['decoy_prefix'])


def run_comet(input_filename, output_psm_filename, database,
              mass_tol_ppm=40,
              output_file_to_fetch='comet.target.txt',
              crux_param_file=None):
    """
    Uses a temporary directory to hold the entire Comet output, from which
    comet.target.txt is copied to output_psm_filename
    """
    if crux_param_file is None or crux_param_file == '':
        param_file = ''
    else:
        param_file = '--parameter-file ' + crux_param_file

    out_dir = Path(tempfile.mkdtemp())
    cmd = (f'crux comet {param_file} --peptide_mass_units 2 '
           f'--peptide_mass_tolerance {mass_tol_ppm} --overwrite T '
           f'--output-dir {out_dir} {input_filename} {database}')

    _run_cmd(cmd, 'Running comet:')
    try:
        shutil.move(str(out_dir / output_file_to_fetch), output_psm_filename)
    except Exception as e:
        print('=' * 20)
        print(output_psm_filename)
        print('=' * 20)
        raise e


def run_percolator(input_targets_filename, input_decoys_filename, output_psm_filename,
                   fdr, decoy_prefix):
    """
    Uses a temporary directory to hold the entire Comet output, from which
    percolator.target.peptides.txt is copied to output_psm_filename
    """
    out_dir = Path(tempfile.mkdtemp())
    cmd = (f'crux percolator --percolator-seed 123 '
           f'--overwrite T --pepxml-output T --mzid-output T '
           f'--output-dir {out_dir} --decoy-prefix {decoy_prefix} --test-fdr {fdr} '
           f'{input_targets_filename} {input_decoys_filename}')

    _run_cmd(cmd, "Running percolator:")
    shutil.move(str(out_dir / 'percolator.target.peptides.txt'), output_psm_filename)


def run_msgf(input_filename, output_filename, config):
    msgf_jar_path = os.getenv('MSGF_JAR_PATH', default='')
    if not msgf_jar_path:
        raise Exception('MSGF_JAR_PATH not set')

    cmd = (f'java -Xmx3500M -jar {msgf_jar_path} '
           f'-s {input_filename} '
           f'-d {config["database"]} '
           f'-tda 1 -decoy {config["decoy_prefix"][:-1]} '
           f'-t {config["mass_tol_ppm"]}ppm '
           f'-inst 2 '
           f'-thread {config["msgf_threads"]} ')

    if config['msgf_modifications']:
        cmd = cmd + f'-mod {config["msgf_modifications"]} '

    _run_cmd(cmd, 'Running MS-GF+:')

    cmd = (f'java -Xmx3500M -cp {msgf_jar_path} edu.ucsd.msjava.ui.MzIDToTsv '
           f'-i {str(Path(input_filename).with_suffix(""))}.mzid -unroll 1')
    _run_cmd(cmd, 'Converting MS-GF+ mzid to tsv:')

    shutil.move(str(Path(input_filename).with_suffix('.tsv')), output_filename)


# TODO: Refactor
def run_msgf_mzid(input_filename, output_filename, config):
    msgf_jar_path = os.getenv('MSGF_JAR_PATH', default='')
    if not msgf_jar_path:
        raise Exception('MSGF_JAR_PATH not set')

    cmd = (f'java -Xmx3500M -jar {msgf_jar_path} '
           f'-s {input_filename} '
           f'-d {config["database"]} '
           f'-tda 1 -decoy {config["decoy_prefix"][:-1]} '
           f'-t {config["mass_tol_ppm"]}ppm '
           f'-inst 2 '
           f'-thread {config["msgf_threads"]} ')

    if config['msgf_modifications']:
        cmd = cmd + f'-mod {config["msgf_modifications"]} '

    _run_cmd(cmd, 'Running MS-GF+:')
    #shutil.move(str(Path(input_filename).with_suffix('.mzid')), output_filename)



def _run_cmd(cmd, log_message):
    logger.info(log_message + ' ' + cmd)
    res = subprocess.run(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logger.info(res.stdout.decode('utf8'))
    logger.info(res.stderr.decode('utf8'))


process_wrapper = {
    'msgf': run_msgf,
    'msgf_mzid': run_msgf_mzid
}


def get_args():
    desc = f'Wrapper module for various MS processes. Available calls:\n{process_wrapper}'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-p',
                        '--process',
                        required=True,
                        type=str,
                        help='Input (process-specific)')
    parser.add_argument('-i',
                        '--input',
                        required=True,
                        type=str,
                        help='Input (process-specific)')
    parser.add_argument('-o',
                        '--output',
                        required=True,
                        type=str,
                        help='Output (process-specific)')
    parser.add_argument('-c',
                        '--config',
                        required=True,
                        type=str,
                        help='YAML experiment config file')
    args = parser.parse_args()
    with open(args.config, 'r') as config_file:
        args.config = yaml.load(config_file, Loader=yaml.FullLoader)
    return args


def main():
    args = get_args()
    process_wrapper[args.process](args.input, args.output, args.config)


if __name__ == '__main__':
    main()
