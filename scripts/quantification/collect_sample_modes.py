import argparse
import inspect
import logging
import os.path as path
import sys

from functools import partial
from multiprocessing import Pool, cpu_count
from pathlib import Path

import coloredlogs
import numpy as np
import pandas as pd
from pyarrow import feather
import scipy.stats as stats
import torch
import yaml


_this_filename = inspect.getframeinfo(inspect.currentframe()).filename
_this_path = Path(_this_filename).parent.resolve()
sys.path.append(str(_this_path.parent))

from util import msproc
from parafac import models

logging.basicConfig(format=msproc.LOG_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)
coloredlogs.install(fmt=msproc.LOG_FORMAT, level='INFO', logger=logger)


def main():
    args = get_args()
    if args.model_set == 'best':
        model_params = get_params_for_best_models(args)
    else:
        model_params = models.read_index(args.config['model_index'],
                                         file_format='feather')

    models_info = [minfo[1] for minfo in list(model_params.iterrows())]

    with Pool(processes=cpu_count()) as pool:
        sample_modes = pool.map(load_model_sample_mode_from_metainfo,
                                models_info)
    sample_modes = [res_df for res_df in sample_modes if res_df is not None]
    logger.info('Read sample modes from best models')

    sample_modes = pd.concat(sample_modes, ignore_index=True)
    feather.write_feather(sample_modes, args.output_file)
    logger.info('Wrote %s' % str(args.output_file))

    save_spectra_and_sample_mode_values(args, sample_modes)


def save_spectra_and_sample_mode_values(args, sample_modes):
    """
    Get the correspondence between decomposed spectra ID (`<scan>` in the mzXML file)
    to the `sample_num` of their corresponding sample (abundance) mode,

    then save it (as a CSV) to config['spectra_with_sample_abundance_file']
    in the experiment directory,
    or 'spectra_with_sample_abundance.csv' if the former is not specified
    """
    spectrum_index = feather.read_feather(args.config['spectrum_index'])
    spectra_with_samples = pd.merge(sample_modes,
                                    spectrum_index,
                                    left_on = ['model_id', 'comp_num'],
                                    right_on = ['model_id', 'spectrum_num'])
    output_fname = args.config.get('spectra_with_sample_abundance_file',
                                   'spectra_with_sample_abundance.csv')

    spectra_with_samples[['scan', 'sample_num', 'abundance']].to_csv(
        args.root_dir / output_fname, index=False
    )
    logger.info('Wrote spectrum-sample abundance values to: ' +
                str(args.root_dir / output_fname))


def load_model_sample_mode_from_metainfo(model_info):
    model_path = model_info['path']
    try:
        model = torch.load(model_path, map_location='cpu')
    except:
        return None
    sample_mode = pd.DataFrame(model[0].numpy(), dtype=np.float32)
    sample_mode = sample_mode.stack().reset_index().rename(index=str,
                                                           columns={'level_0': 'sample_num',
                                                                    'level_1': 'comp_num',
                                                                     0: 'abundance'})
    sample_mode = sample_mode.assign(model_id = model_info['model_id'])

    abundance_cv_per_comp = sample_mode.groupby('comp_num')['abundance'].apply(stats.variation)
    abundance_cv_per_comp.name = 'cv_sample_mode'
    sample_mode =  pd.merge(sample_mode, abundance_cv_per_comp,
                            left_on='comp_num', right_on=abundance_cv_per_comp.index)

    return sample_mode


def get_params_for_best_models(args):
    best_models_fname = args.root_dir / args.config['best_models']
    if not path.exists(best_models_fname):
        logger.error('Stopping: Could not find file listing best PARAFAC models: ' +
                     str(best_models_fname))
        sys.exit(1)
    model_params = pd.read_csv(best_models_fname)
    model_params['path'] = model_params.apply(partial(_model_path_from_params,
                                                      args.config['slices_location']),
                                              axis='columns')
    model_index = models.read_index(args.config['model_index'], file_format='feather')
    model_params = pd.merge(model_index['model_id'], model_params,
                            how='right', on=['model_id'])
    return model_params


def _model_path_from_params(slies_location, model_params):
    model_path = path.join(slies_location,
                           f"swath_lower_adjusted={model_params['swath_start']:.2f}",
                           f"rt_window={model_params['rt_window']:.1f}",
                           f"parafac_model_F{int(model_params['ncomp'])}.pt")
    return model_path


def get_args():
    desc = 'Quantify peptides from best PARAFAC models'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-c',
                        '--config',
                        required=True,
                        type=str,
                        help='YAML experiment config file')
    parser.add_argument('--model_set',
                        required=False,
                        type=str,
                        default='best',
                        help='Set of models to collect sample modes from. '
                             'Values: best, all  Default: %(default)s)')
    args = parser.parse_args()
    with open(args.config, 'r') as config_file:
        args.config = yaml.load(config_file, Loader=yaml.FullLoader)
    if args.config['analysis_pipeline'] not in ['crux', 'tpp', 'msgf+']:
        logger.error("Analytic pipeline " + args.config['analysis_pipeline'] + ' not supported')
        raise Exception('Invalid pipeline config value: analysis_pipeline = ' + args.config['analysis_pipeline'])
    args.root_dir = Path(args.config['root_dir'])
    args.output_file = args.root_dir / f'sample_modes_{args.model_set}_models.feather'
    return args


if __name__ == '__main__':
    main()
