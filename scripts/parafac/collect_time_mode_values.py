import argparse
import inspect
import logging
import os
import os.path as path
import sys

from functools import partial
from multiprocessing import Pool, cpu_count
from pathlib import Path

import coloredlogs
import numpy as np
import pandas as pd
import pyarrow.feather as feather
import scipy.signal as sig
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
    model_params = models.read_index(Path(args.config['root_dir']) / args.config['model_index'],
                                     file_format='feather')
    models_info = [minfo[1] for minfo in list(model_params.iterrows())]

    with Pool(processes=cpu_count()) as pool:
        model_peak_count = pool.map(partial(get_model_time_mode_peak_counts, args),
                                    models_info)
    model_peak_count = [res_df for res_df in model_peak_count if res_df is not None]
    logger.info('Read time modes from all models')

    model_peak_count = pd.concat(model_peak_count, ignore_index=True)

    model_peak_count_filename = Path(args.config['root_dir']) / args.config['time_modes_values']
    if not model_peak_count_filename.parent.exists():
        os.makedirs(model_peak_count_filename.parent)
    feather.write_feather(model_peak_count, model_peak_count_filename)
    logger.info(f"Wrote {Path(args.config['root_dir']) / args.config['time_modes_values']}")


def get_model_time_mode_peak_counts(args: argparse.Namespace,
                                    model_info: pd.DataFrame) -> pd.DataFrame:
    model_path = model_info['path']
    try:
        model = torch.load(model_path, map_location='cpu')
    except:
        logger.warning('Could not load model ' + model_path)
        return None

    time_mode = model[1].numpy()
    avg_peak_window_frac = (args.config['avg_peak_fwhm_sec']
                            / args.config['window_size_sec'])
    expected_peak_width = time_mode.shape[0] * avg_peak_window_frac

    npeaks = np.apply_along_axis(partial(count_peaks, expected_peak_width),
                                 axis=0,
                                 arr=time_mode)

    npeaks = pd.DataFrame.from_records(enumerate(npeaks), columns=['comp_num', 'npeaks'])
    npeaks['comp_num'] = npeaks['comp_num'].astype(np.uint8)
    npeaks['npeaks'] = npeaks['npeaks'].astype(np.uint8)
    npeaks = npeaks.assign(model_id=model_info['model_id'])
    return npeaks


def count_peaks(expected_peak_width: float, profile: np.array) -> int:
    # Remove very low background values to help peak detection
    clipped_profile = np.copy(profile)
    null_threshold = clipped_profile.max() * 0.1
    clipped_profile[clipped_profile <= null_threshold] = 0

    peaks = sig.find_peaks_cwt(clipped_profile,
                               widths=np.arange(1, expected_peak_width * 2))
    return peaks.shape[0]


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
    args = parser.parse_args()
    with open(args.config, 'r') as config_file:
        args.config = yaml.load(config_file, Loader=yaml.FullLoader)
    return args


if __name__ == '__main__':
    main()
