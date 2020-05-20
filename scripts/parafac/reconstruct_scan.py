"""
Reconstruct m/z x RT matrix from PARAFAC ouput, for a single scan
"""
import argparse
import inspect
import json
import logging
import os.path as path
import sys

from functools import partial
from multiprocessing import Pool, cpu_count
from pathlib import Path

import coloredlogs
import numpy as np
from numpy import unravel_index
import pandas as pd
import pyarrow.feather as feather
import torch
import yaml


_this_filename = inspect.getframeinfo(inspect.currentframe()).filename
_this_path = Path(_this_filename).parent.resolve()
sys.path.append(str(_this_path.parent))

from ..util import msproc
import models

logging.basicConfig(format=msproc.LOG_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)
coloredlogs.install(fmt=msproc.LOG_FORMAT, level='INFO', logger=logger)


def main():
    args = get_args()
    best_model_params = models.read_index(args.config['best_models']).reset_index()

    with Pool(processes=cpu_count()) as pool:
        slice_maps = pool.map(partial(reconstruct_slice_map, args),
                                      best_model_params[['path', 'model_id']].values)

    slice_maps = [reconstr_map for reconstr_map in slice_maps if reconstr_map is not None]
    logger.info('Reconstructed maps for all slices')

    slice_maps = pd.concat(slice_maps, ignore_index=True)
    feather.write_feather(slice_maps, args.output)
    logger.info(f'Wrote {args.output}')


def reconstruct_slice_map(args: argparse.Namespace, model_info: list) -> pd.DataFrame:
    """
    Sum component traces for a single scan, i.e.
    sum( ( mass[k] x time[k] ) * sample[k][scan] )
    """
    model_path = model_info[0]
    model_id = model_info[1]
    try:
        model = torch.load(model_path, map_location='cpu')
        sample_mode = model[0].numpy()
        time_mode = model[1].numpy()
        mass_mode = model[2].numpy()
    except:
        logger.warning('Could not load model ' + model_path)
        return None

    # Load tensor information
    prop_file_path = Path(model_path).parent / 'tensor_properties.json'
    with open(prop_file_path, 'r') as slice_tensor_properties_file:
        slice_tensor_properties = json.load(slice_tensor_properties_file)
        mz_indices = msproc.expand_list_of_mz_indices_to_df(slice_tensor_properties['mz_indices'])
        mz_indices = mz_indices.assign(mz_idx=mz_indices.index)
        sample_fnames = slice_tensor_properties['samples']

    # Fix for shuffled filenames in existing results
    # No effect on properly sorted filenames in slice
    scan_fname = sorted(sample_fnames)[args.scan_num]
    scan_idx_in_tensor = sample_fnames.index(scan_fname)

    # Select only unimodal components
    npeaks = feather.read_feather(Path(args.config['best_models']).parent / 'time_mode_values_all_models.feather')
    component_numbers = npeaks[(npeaks['model_id'] == model_id) & (npeaks['npeaks'] == 1)]['comp_num'].values
    del npeaks  # try to free some memory

    if component_numbers.shape[0] == 0:
        logger.info('Slice model had no unimodal components: ' + model_path)
        return None

    # Cross product and sample scaling for all components
    reconstr_map = np.zeros((time_mode.shape[0], mass_mode.shape[0]))

    for k in component_numbers:
        trace = np.outer(time_mode[:, k], mass_mode[:, k])
        trace *= sample_mode[scan_idx_in_tensor, k]
        reconstr_map = np.add(reconstr_map, trace)

    # Scaling: Get intensity from slice matching the reconstructed map's max value
    slice_tensor = np.load(Path(model_path).parent / 'slice_tensor.npy')
    argmax_reconstr_x, argmax_reconstr_y = unravel_index(reconstr_map.argmax(), reconstr_map.shape)
    if argmax_reconstr_x >= slice_tensor[scan_idx_in_tensor].shape[0]:
        argmax_reconstr_x = slice_tensor[scan_idx_in_tensor].shape[0] - 1
        logger.warning('Argmax out of bounds: clipping')
    if argmax_reconstr_y >= slice_tensor[scan_idx_in_tensor].shape[1]:
        argmax_reconstr_y = slice_tensor[scan_idx_in_tensor].shape[1] - 1
        logger.warning('Argmax out of bounds: clipping')
    slice_value_at_argmax = slice_tensor[scan_idx_in_tensor][(argmax_reconstr_x, argmax_reconstr_y)]

    # Scaling: Get model R^2 to scale reconstructed map
    rsq_vals = pd.read_csv(args.config['rsq_best_models'])
    try:
        model_rsq = rsq_vals[rsq_vals['model_id'] == model_id]['Rsq'].values[0]
    except IndexError:
        model_rsq = 0.5
        logger.warning('R^2 value missing for model. Using R^2 = 0.5')

    scale_ratio = (slice_value_at_argmax * model_rsq) / reconstr_map.max()
    reconstr_map *= scale_ratio

    # Convert to long format
    reconstr_map = pd.DataFrame(reconstr_map.transpose()).reset_index()
    reconstr_map = reconstr_map.melt(id_vars='index', value_vars=None,
                                     var_name='cycle', value_name='intensity')
    reconstr_map = pd.merge(reconstr_map, mz_indices, left_on='index', right_on='mz_idx')
    reconstr_map = reconstr_map.drop(columns=['index', 'mz_idx'])
    reconstr_map = reconstr_map.assign(model_id=model_id)

    # Reduce size
    reconstr_map['cycle'] = reconstr_map['cycle'].astype('uint16')
    reconstr_map['model_id'] = reconstr_map['model_id'].astype('uint32')
    reconstr_map['intensity'] = reconstr_map['intensity'].astype('float32')

    logger.info('Reconstructed ' + model_path)
    return reconstr_map


def get_args():
    desc = 'Reconstruct m/z x RT matrix from PARAFAC ouput, for a single scan'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-c',
                        '--config',
                        required=True,
                        type=str,
                        help='YAML experiment config file')
    parser.add_argument('-o',
                        '--output',
                        required=True,
                        type=str,
                        help='Output file')
    parser.add_argument('-s',
                        '--scan_num',
                        required=False,
                        type=int,
                        default=0,
                        help='Scan number. Default:  %(default)s)')
    args = parser.parse_args()
    with open(args.config, 'r') as config_file:
        args.config = yaml.load(config_file, Loader=yaml.FullLoader)
    return args


if __name__ == '__main__':
    main()
