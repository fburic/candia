import argparse
import inspect
import locale
import logging
import os
from functools import partial
from pathlib import Path

import parse
import sys
from multiprocessing import Pool, cpu_count

import coloredlogs
import numpy as np
import pandas as pd
import pyarrow.feather as feather
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
    slice_paths = get_slice_paths_from_model_index(args)
    if not args.skip_id_results:
        save_concatenated_id_results_from_all_slices(slice_paths, args)
    if not args.skip_decomp_measures:
        save_decomposition_measures_from_all_slices(slice_paths, args)


def get_slice_paths_from_model_index(args):
    model_index = models.read_index(args.config['model_index'], file_format='feather')
    slice_paths = model_index['path'].values
    slice_paths = map(lambda path: os.path.dirname(path), slice_paths)
    slice_paths = list(set(slice_paths))
    return slice_paths


def save_concatenated_id_results_from_all_slices(slice_paths, args):
    logger.info('Concatenating identification CSV results from all slices...')
    result_filepaths = map(lambda path: os.path.join(path, 'id_results.csv'),
                           slice_paths)
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(partial(read_csv_result_file_as_dataframe,
                                   add_slice_params=False,
                                   hide_warnings=args.hide_warnings),
                           result_filepaths)

    results_df = [res_df for res_df in results if res_df is not None]
    if results_df:
        logger.info('... Read %d result files ...' % len(results_df))
        if args.config['analysis_pipeline'] == 'crux':
            collected_results_file = args.config['comet_target_results_file']
        elif args.config['analysis_pipeline'] == 'tpp':
            collected_results_file = args.config['tadnem_results_file']

        results = pd.concat(results_df)
        results.to_csv(collected_results_file, index=False, compression='gzip')
        logger.info('... done. Wrote ' + collected_results_file)
    else:
        logger.error('... There is no ID output data in any slice!')


def save_decomposition_measures_from_all_slices(slice_paths, args):
    logger.info('Concatenating decomposition measures from all slices...')
    decomp_measures_filepaths = map(lambda path: os.path.join(path, 'decomp_measures.csv'),
                                    slice_paths)
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(partial(read_csv_result_file_as_dataframe,
                                   add_slice_params=True,
                                   hide_warnings=args.hide_warnings),
                           decomp_measures_filepaths)

    results_df = [res_df for res_df in results if res_df is not None]
    if results_df:
        logger.info('... Read %d result files ...' % len(results_df))
        results = pd.concat(results_df)
        feather.write_feather(results, args.config['decomp_measures'])
        logger.info('... done. Wrote ' + args.config['decomp_measures'])
    else:
        logger.error('... No decomposition measures found!')


def read_csv_result_file_as_dataframe(filename, add_slice_params=False, hide_warnings=False):
    try:
        result_table = pd.read_csv(filename)

        if add_slice_params:
            swath_start, rt_window = extract_slice_params_from_path(filename)
            if swath_start is not None:
                result_table = result_table.assign(swath_start=swath_start,
                                                   rt_window=rt_window)

        result_table['swath_start'] = result_table['swath_start'].apply(
            lambda s: format(s, '.2f')
        )
        result_table['rt_window'] = result_table['rt_window'].astype(np.uint8)
        result_table['ncomp'] = result_table['ncomp'].astype(np.uint8)
        return result_table

    except pd.errors.EmptyDataError as e:
        if not hide_warnings:
            logger.warning('Skipping empty ID output file: ' + filename)
            logger.warning(sys.exc_info()[0])
            logger.warning(e)
        return None
    except FileNotFoundError as e:
        if not hide_warnings:
            logger.warning('Skipping missing slice ID output file: ' + filename)
            logger.warning(sys.exc_info()[0])
            logger.warning(e)
        return None


def extract_slice_params_from_path(filepath):
    """
    Parse a path like .../swath_lower_adjusted=500.5/rt_window=12.0/file.txt
    and return slice parameters (swath_start, rt_window) as numerical values.
    """
    swath_start = None
    rt_window = None
    locale.setlocale(locale.LC_NUMERIC, 'en_US.UTF-8')
    values = list(parse.findall("={}/", filepath))
    if values:
        swath_start = locale.atof(list(values[0])[0])
        rt_window = locale.atof(list(values[1])[0])
    return swath_start, rt_window


def get_args():
    desc = 'Collect model evaluation results from all slices'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-c',
                        '--config',
                        required=True,
                        type=str,
                        help='YAML experiment config file')
    parser.add_argument('-o',
                        '--pep_counts_file',
                        required=False,
                        type=str,
                        help='File in which to save high-confidence peptide counts')
    parser.add_argument('--hide-warnings',
                        required=False,
                        action='store_true',
                        help='Hide warnings due to missing slice data')
    parser.add_argument('--skip-id-results',
                        required=False,
                        action='store_true',
                        help='')
    parser.add_argument('--skip-decomp-measures',
                        required=False,
                        action='store_true',
                        help='')

    args = parser.parse_args()
    with open(args.config, 'r') as config_file:
        args.config = yaml.load(config_file)

    if args.config['analysis_pipeline'] not in ['crux', 'tpp']:
        logger.error("Analytic pipeline " + args.config['analysis_pipeline'] + ' not supported')
        raise Exception('Invalid pipeline config value: analysis_pipeline = ' + args.config['analysis_pipeline'])

    return args


if __name__ == '__main__':
    main()