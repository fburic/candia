"""
Pick most promising models as assessed by Comet,
concatenate into single mzXML file,
and use Comet + Percolator to get high-confidence PSMs
"""
import argparse
import difflib
import os
from collections import namedtuple
from functools import partial
import inspect
import json
import logging
import os.path as path
import subprocess
import sys
from pathlib import Path

import coloredlogs
import numpy as np
import pandas as pd
import torch
import yaml
from tqdm import tqdm

_this_filename = inspect.getframeinfo(inspect.currentframe()).filename
_this_path = Path(_this_filename).parent.resolve()
sys.path.append(str(_this_path.parent))

from ..util import msproc
from ..parafac import models

logging.basicConfig(format=msproc.LOG_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)
coloredlogs.install(fmt=msproc.LOG_FORMAT, level='INFO', logger=logger)


def main():
    args = get_args()

    if not args.skip_xml_generation:
        concatenate_best_models_to_mzxml(args)

    if args.config['analysis_pipeline'] == 'crux':
        comet_target_dir = path.join(args.config['best_models_crux_out_dir'], 'comet_target')
        comet_decoy_dir = path.join(args.config['best_models_crux_out_dir'], 'comet_decoy')

        get_comet_psm(args, comet_target_dir, args.config['database'])
        get_comet_psm(args, comet_decoy_dir, args.config['decoy_database'])
        get_percolator_separation(args, comet_target_dir, comet_decoy_dir)

    elif args.config['analysis_pipeline'] == 'msgf+':
        get_msgf_psm(args)

    else:
        raise NotImplementedError


def concatenate_best_models_to_mzxml(args):
    """
    Convert mass modes of best models into <scan> entries and write to a mzXML file
    """
    logger.info('Concatenating best models to single mzXML file with unique scan IDs...')
    model_params = get_params_for_best_models(args).to_dict('records')

    spectrum_index = models.read_index(args.config['spectrum_index'], file_format='feather')
    spectrum_index = spectrum_index.drop(columns='model_id')  # Save some memory
    isol_windows = pd.read_csv(args.config['swath_windows_adjusted'])
    intensity_cutoff_bin = args.config['intensity_lower_percentage_cutoff']

    with msproc.MassMode2MzXMLEncoder(args.config['best_models_mzxml']) as mzxml_encoder:

        for model_info in tqdm(model_params):
            try:
                prop_file_path = path.join(path.dirname(model_info['path']), 'tensor_properties.json')
                with open(prop_file_path, 'r') as slice_tensor_properties_file:
                    slice_tensor_properties = json.load(slice_tensor_properties_file)
                    try:
                        mz_indices = slice_tensor_properties['mz_indices']
                    except KeyError:
                        logger.warning('Could not properly read tensor properties %s, mz_indices not found. Skipping.'
                                       % prop_file_path)
                        continue

                window = isol_windows[np.isclose(isol_windows['swath_lower_adjusted'],
                                                 model_info['swath_start'])]
                isolation_window_center = \
                    (window['swath_lower_adjusted'] + window['swath_upper_adjusted']) / 2
                isolation_window_center = isolation_window_center.values[0]

                slice_model = torch.load(model_info['path'], map_location='cpu')
                mass_mode = slice_model[2].numpy()
                slice_model = None  # Save some memory

                # def mass_mode_ndarray_as_df_indexed_by_mz():
                mass_mode = pd.DataFrame(mass_mode, dtype=np.float32)
                mass_mode = mass_mode.assign(mz_indices=mz_indices)

                # Replace ordinal spectrum (component) number with globally unique scan ID
                spectrum_to_scan_for_model = models.spectrum_index_entries_for_model(model_info, spectrum_index)
                spectrum_to_scan_for_model = spectrum_to_scan_for_model[['spectrum_num', 'scan']].values
                spectrum_to_scan_for_model = {spectrum_to_scan[0]: spectrum_to_scan[1]
                                              for spectrum_to_scan in spectrum_to_scan_for_model}
                mass_mode = mass_mode.rename(index=str, columns=spectrum_to_scan_for_model)

                mass_mode = mass_mode.set_index('mz_indices')
                mass_mode = msproc.as_dataframe_with_mz_index(mass_mode,
                                                              mass_mode.index.values)

                mzxml_encoder.write_scans_from_mass_mode_df(mass_mode,
                                                            intensity_cutoff_bin=intensity_cutoff_bin,
                                                            isolation_window_center=isolation_window_center)
                logger.info('Wrote scans from model ' + model_info['path'])

            except FileNotFoundError as e:
                logger.warning(sys.exc_info()[0])
                logger.warning(e)

    logger.info('... Export to mzXML done.')


def get_params_for_best_models(args):
    model_params = pd.read_csv(args.config['best_models'])
    model_params['path'] = model_params.apply(partial(_model_path_from_params,
                                                      args.config['slices_location']),
                                              axis='columns')
    return model_params


def _model_path_from_params(slies_location, model_params):
    model_path = path.join(slies_location,
                           f"swath_lower_adjusted={model_params['swath_start']:.2f}",
                           f"rt_window={model_params['rt_window']:.1f}",
                           f"parafac_model_F{int(model_params['ncomp'])}.pt")
    return model_path


def get_comet_psm(args, comet_output_dir, database):
    run_comet(output_dir=comet_output_dir,
              mass_mode_filename=args.config['best_models_mzxml'],
              library=database,
              mass_tol_ppm=args.config['mass_tol_ppm'],
              crux_param_file=None)


def get_percolator_separation(args, comet_target_dir, comet_decoy_dir):
    comet_target_results = comet_target_dir + '/comet.target.txt'
    comet_decoy_results = comet_decoy_dir + '/comet.target.txt'
    percolator_out_dir = args.config['best_models_crux_out_dir']
    run_percolator(percolator_out_dir,
                   args.config['percolator_fdr'],
                   comet_target_results,
                   comet_decoy_results,
                   decoy_prefix=args.config['decoy_prefix'])


def run_comet(output_dir, mass_mode_filename, library, mass_tol_ppm=40, crux_param_file=None):
    # TODO: Move this to some IO / util package and refactor uses in other scripts
    if crux_param_file is None:
        crux_param_file = ''
    else:
        crux_param_file = '--parameter-file ' + crux_param_file

    cmd = 'crux comet {param_file} --peptide_mass_units 2 --peptide_mass_tolerance {tol} --overwrite T --output-dir {out_dir} {in_file} {lib}'
    cmd = cmd.format(param_file=crux_param_file,
                     tol=mass_tol_ppm,
                     out_dir=output_dir,
                     in_file=mass_mode_filename,
                     lib=library)
    _run_cmd(cmd, 'Running comet:')


def run_percolator(output_dir, fdr, comet_target_results, comet_decoy_results, decoy_prefix='reverse_'):
    # TODO: Move this to some IO / util package and refactor uses in other scripts
    cmd = ('crux percolator --percolator-seed 123 '
           '--overwrite T --pepxml-output T --mzid-output T '
           '--output-dir {out_dir} --decoy-prefix {decoy_prefix} --test-fdr {fdr} '
           '{comet_target_results} {comet_decoy_results}')
    cmd = cmd.format(out_dir=output_dir,
                     decoy_prefix=decoy_prefix,
                     fdr=fdr,
                     comet_target_results=comet_target_results,
                     comet_decoy_results=comet_decoy_results)
    _run_cmd(cmd, "Running percolator:")


def get_msgf_psm(args):
    msgf_jar_path = os.path.expanduser(os.getenv('MSGF_JAR_PATH', default=''))
    if not msgf_jar_path:
        raise Exception('MSGF_JAR_PATH not set')

    best_models_path = Path(args.config["best_models_mzxml"])
    cmd = (f'java -Xmx3500M -jar {msgf_jar_path} '
           f'-s {best_models_path} '
           f'-d {args.config["database"]} '
           f'-tda 1 -decoy {args.config["decoy_prefix"][:-1]} '
           f'-t {args.config["mass_tol_ppm"]}ppm '
           f'-inst 2 '
           f'-mod {args.config["msgf_modifications"]}')
    _run_cmd(cmd, 'Running MS-GF+:')

    cmd = (f'java -Xmx3500M -cp {msgf_jar_path} edu.ucsd.msjava.ui.MzIDToTsv '
           f'-i {best_models_path.with_suffix("")}.mzid -unroll 1')
    _run_cmd(cmd, 'Converting MS-GF+ mzid to tsv:')


def _run_cmd(cmd, log_message):
    logger.info(log_message + ' ' + cmd)
    res = subprocess.run(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logger.info(res.stdout.decode('utf8'))
    logger.info(res.stderr.decode('utf8'))


def get_args():
    desc = 'Identify spectra in PARAFAC models using comet and percolator'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-c',
                        '--config',
                        required=True,
                        type=str,
                        help='YAML experiment config file')
    parser.add_argument('--skip_xml_generation',
                        required=False,
                        action='store_true',
                        help='Skip creating mzXML file from best models')
    args = parser.parse_args()
    with open(args.config, 'r') as config_file:
        args.config = yaml.load(config_file, Loader=yaml.FullLoader)
    if args.config['analysis_pipeline'] not in ['crux', 'tpp', 'msgf+']:
        logger.error("Analytic pipeline " + args.config['analysis_pipeline'] + ' not supported')
        raise Exception('Invalid pipeline config value: analysis_pipeline = ' + args.config['analysis_pipeline'])
    return args


def get_test_args():
    Args = namedtuple('Args', ['config'])
    args = Args(config={"slices_location": "test/models/slices",
                        "analysis_pipeline": "crux",
                        "spectrum_index": "test/models/spectrum_index.h5",
                        "swath_windows_adjusted": "test/models/swaths_adjusted.csv",
                        "best_models": "test/models/best_models_acc_to_comet.csv",
                        "best_models_mzxml": "/tmp/best_models_acc_to_comet.mzXML",
                        "best_models_crux_out_dir": "test/models/crux-output-concat",
                        "percolator_fdr": 0.01})
    return args


def test_concatenate_best_models_to_mzxml():
    # FIXME: mass mode is too large. M/z indices are a small list.

    expected_mzxml = 'test/models/best_models_acc_to_comet.mzXML'
    args = get_test_args()
    concatenate_best_models_to_mzxml(args)

    with open(args.config['best_models_mzxml'], 'r') as test:
        with open(expected_mzxml, 'r') as expected:
            diff = list(difflib.unified_diff(test.readlines(), expected.readlines(), n=0))

    assert len(diff) == 0


if __name__ == '__main__':
    main()
