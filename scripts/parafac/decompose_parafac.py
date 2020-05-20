import argparse
import inspect
import json
import logging
import os
import sys
from pathlib import Path

import coloredlogs
import numpy as np
import pandas as pd

from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve

import torch
import torch.multiprocessing as mp
import tensorly as tl
from tensorly.decomposition.candecomp_parafac import initialize_factors
from tensorly import backend as T
from tensorly.base import unfold
from tensorly.kruskal_tensor import kruskal_to_tensor
from tensorly.tenalg import khatri_rao

_this_filename = inspect.getframeinfo(inspect.currentframe()).filename
_this_path = Path(_this_filename).parent.resolve()
sys.path.append(str(_this_path.parent))

LOG_FORMAT = "[%(asctime)s] [PID %(process)d] %(levelname)s:\t%(filename)s:%(funcName)s():%(lineno)s:\t%(message)s"

logging.basicConfig(format=LOG_FORMAT)
logger = logging.getLogger(__name__)
coloredlogs.install(fmt=LOG_FORMAT, level='INFO', logger=logger)


def init_frameworks(args):
    # These are needed to talk with the CUDA MPS server
    # https://docs.nvidia.com/deploy/pdf/CUDA_Multi_Process_Service_Overview.pdf
    os.environ['pipe_id'] = str(args.pipe_id)
    os.environ['CUDA_MPS_PIPE_DIRECTORY'] = os.path.expandvars('${TMPDIR}/mps-${pipe_id}/pipe')
    os.environ['CUDA_MPS_LOG_DIRECTORY'] = os.path.expandvars('${TMPDIR}/mps-${pipe_id}/log')
    os.unsetenv('CUDA_VISIBLE_DEVICES')

    tl.set_backend(args.backend)
    override_tensorly_methods(args)
    set_random_seed(args.rnd_seed, args)


def main():
    # Only way this will run. So no error handling: die if this doesn't work.
    mp.set_start_method('spawn')

    args = get_args()
    model_filenames = {ncomp: os.path.join(args.output_dir,
                                           'parafac_model_F{}.pt'.format(ncomp))
                       for ncomp in args.component_range}

    try:
        tensor = load_tensor(args.input_tensor_file)
    except TrivialTensorlException as e:
        logger.warning(str(e))
        with (Path(args.output_dir) / 'decomp_measures.csv').open(mode='w') as logfile:
            logfile.write('ncomp,Rsq,iterations\n')
        return 0

    tensor = impute_missing_values(tensor)
    init_frameworks(args)
    tensor = to_tensor(tensor, args)

    logger.info('Decomposition range: %d - %d components' % (args.min_components,
                                                             args.max_components))
    models = {}
    for ncomp in args.component_range:
        if args.skip_completed and os.path.exists(model_filenames[ncomp]):
            logger.info('Skipping fpr F = %d Model file already exists.' % ncomp)
            continue

        else:
            models[ncomp] = compute_model(ncomp, tensor, args)

    measures = []
    for ncomp in args.component_range:
        if ncomp in models:
            model, rsq, n_iters = models[ncomp]
            save_model(model, model_filenames[ncomp])
            measures.append((ncomp, rsq, n_iters))

    if args.use_gpu:
        torch.cuda.empty_cache()

    write_decomposition_measures_to_logfile(measures, args, logfile='parafac.log')
    write_decomposition_measures_to_logfile(measures, args, logfile='decomp_measures.csv') # TODO: quick hack


def compute_model(n_components, tensor, args):
    """
    Decompose the tensor for F = n_components and store the model in a pt file.
    """
    results = modif_non_negative_parafac(tensor,
                                         rank=n_components,
                                         n_iter_max=args.max_iter,
                                         init=args.init,
                                         tol=args.tolerance,
                                         random_state=args.rnd_seed,
                                         verbose=False)
    model, reconcstr_error, n_iters = results
    rsq = (1 - reconcstr_error ** 2)

    logger.info('Decomposed for F = %d\t: R^2 = %f\titer = %d'
                % (n_components, rsq, n_iters))
    return (model, rsq, n_iters)


def override_tensorly_methods(args):

    if args.backend == 'pytorch':
        def tensor_on_device(data,
                             dtype=tl.float32,
                             device=args.device_id,   # FIX: was 'cpu'
                             requires_grad=False):
            """Override Tensor creation to use our device (ver 0.4.3)"""
            if isinstance(data, np.ndarray):
                return torch.tensor(data.copy(), dtype=dtype, device=device, requires_grad=requires_grad)
            return torch.tensor(data, dtype=dtype, device=device, requires_grad=requires_grad)

        tl.tensor = tensor_on_device


def impute_missing_values(tensor):
    logger.info('Imputing missing values ...')
    imputed_slices = []
    for chromatogram in tensor:
        imputed_slices.append(impute_missing_values_for_chromatogram(chromatogram))
    logger.info('... done.')
    return np.stack(imputed_slices)


def impute_missing_values_for_chromatogram(chromatogram):
    gaussian_kernel = Gaussian1DKernel(stddev=0.5, x_size=5)
    imputed_map = np.apply_along_axis(blur_chromatorgram, 0, chromatogram, (gaussian_kernel))
    original_vals_idx = np.isfinite(chromatogram)
    imputed_map[original_vals_idx] = chromatogram[original_vals_idx]
    imputed_map = np.nan_to_num(imputed_map, copy=False)
    imputed_map = np.clip(imputed_map, a_min=0, a_max=None)
    return imputed_map


def blur_chromatorgram(chromatogram, *args):
    kernel = args[0]
    blurred_chromatogram = convolve(chromatogram, kernel, nan_treatment='fill')
    missing_data_idx = np.isnan(chromatogram)
    chromatogram[missing_data_idx] = blurred_chromatogram[missing_data_idx]
    return chromatogram


def remove_ms1_values(tensor, args):
    with open(args.tensor_properties, 'r') as slice_tensor_properties_file:
        slice_tensor_properties = json.load(slice_tensor_properties_file)

    ms2_indices = mz_indices_for_level(slice_tensor_properties['mz_indices'],
                                       level=2)
    return tensor[:, :, ms2_indices]


def mz_indices_for_level(mz_index_list, level):
    """
    Return list indices for the mz values at the given MS level
    DUPLICATED from msproc for performance
    """
    lvl_marker = '_ms' + str(level)
    return [i for i, mz in enumerate(mz_index_list) if mz.endswith(lvl_marker)]


def convert_tensor_to_df(tensor, mz_indices,
                         sample_axis=0, time_axis=1, mz_axis=2):
    """
    Return a pandas.DataFrame with columns
    ['sample_no', 'cycle', 'mz', 'level', 'intensity']
    """
    # Tensor is processed as (sample x time x m/z)
    if sample_axis != 0 or time_axis != 1 or mz_axis != 2:
        tensor = tensor.transpose(axes=[sample_axis, time_axis, mz_axis])
    sample_slices = np.split(tensor, tensor.shape[0], axis=0)
    mz_col_idx_to_name_mapping = {i: mz for i, mz in enumerate(mz_indices)}
    slice_df_per_sample = []
    for sample_no, sample_slice in enumerate(sample_slices):
        slice_df = pd.DataFrame(sample_slice[0])
        slice_df = slice_df.rename(columns=mz_col_idx_to_name_mapping)
        slice_df['cycle'] = slice_df.index
        slice_df = slice_df.melt(id_vars='cycle', var_name='mz',
                                 value_name='intensity')
        # TODO: Spliting is slow
        mz_level = slice_df['mz'].str.split('_ms', expand=True)
        slice_df['level'] = mz_level[1].astype(np.uint8)
        slice_df['mz'] = mz_level[0]
        slice_df['sample_no'] = sample_no
        slice_df_per_sample.append(slice_df)
    slice_df_per_sample = pd.concat(slice_df_per_sample, ignore_index=True)
    return slice_df_per_sample


def preprocess(tensor, ms1_indices):
    """
    Scale (sample, time) 2D slabs across m/z axis.
    Keep track of maxima in tensor properties file.

    Min-Max scaling deemed appropriate:
    - preserves shape of original distribution (outliers still important)
    - expect a true baseline of zero (no signal),
      so subtracting the min should not remove relevant peaks

    However, at this point we've already filtered out intensities <= 0
    and those coordinates likely to be imputed to 0,
    so this essentially becomes scaling to max.
    """
    # Tensor shape = (sample x rt x m/z)

    scaling_weights = []
    for mz in range(tensor.shape[2]):
        weight = np.sqrt(np.mean(np.square(tensor[:, :, mz])))
        tensor[:, :, mz] /= weight
        scaling_weights.append(weight)

    scaling_weights = np.array(scaling_weights).reshape((len(scaling_weights), 1))
    logger.info('Preprocessed tensor')
    return tensor, scaling_weights


def postprocess(model, scaling_weights, ms1_indices):
    """Scale back mass mode"""
    model[2] = torch.mul(model[2].t(), scaling_weights).t()
    model[2][ms1_indices, :] *= 1000
    return model


def expand_list_of_mz_indices_to_df(mz_indices):
    """
    Take list of string mz_indices as e.g. ['1000.75_ms1', '100.5_ms2']
    and converts to pandas.DataFrame with columns ['mz', 'level']

    DUPLICATED from msproc for performance
    """
    mz_indices = pd.DataFrame(mz_indices)
    mz_indices = mz_indices.assign(mz=mz_indices[0].str.split('_', expand=True)[0].astype(np.float32))
    mz_indices = mz_indices.assign(level=mz_indices[0].str.split('ms', expand=True)[1].astype(np.uint8))
    mz_indices = mz_indices.drop(columns=0)
    return mz_indices


def save_tensor_properties(tensor_properties, slice_dir):
    with open(os.path.join(slice_dir, 'tensor_properties.json'), 'w') as properties_file:
        json.dump(tensor_properties, properties_file, indent=4)


def save_model(model, model_filename):
    torch.save(model, model_filename)
    logger.info('Saved modes to %s' % model_filename)


def get_device_identifier(args):
    if args.backend == 'pytorch':
        if args.use_gpu:
            return 'cuda'
        else:
            return 'cpu'
    else:
        return 'cpu'


def to_tensor(tensor, args):
    if args.backend in ['tensorflow', 'pytorch']:
        return tl.tensor(tensor, device=args.device_id)
    else:
        return tensor


def load_tensor(input_tensor):
    tensor = np.load(input_tensor)
    if trivial_tensor(tensor):
        raise TrivialTensorlException('Tensor contains too little data to be useful.')
    else:
        logger.info('Loaded (sample x rt x m/z) = %s tensor.' % (tensor.shape, ))
        return tensor

class TrivialTensorlException(Exception):
    pass


def trivial_tensor(tensor):
    """These tensors contain too little data to be useful."""
    dim_samples, dim_time, dim_mass = tensor.shape
    if dim_samples < 2 or dim_time < 3 or dim_mass < 3:
        return True
    else:
        return False


def set_random_seed(rnd_seed, args):
    np.random.seed(rnd_seed)
    if args.backend == 'pytorch':
        torch.manual_seed(rnd_seed)
        if args.use_gpu:
            torch.cuda.manual_seed(rnd_seed)
            torch.cuda.manual_seed_all(rnd_seed)


def model_rsq(input_tensor, model_tensor):
    diff = T.norm(input_tensor - kruskal_to_tensor(model_tensor), 2)
    reconcstr_error = diff / T.norm(input_tensor, 2)
    rsq = 1 - reconcstr_error ** 2
    return rsq


def modif_non_negative_parafac(tensor, rank, n_iter_max=100, init='svd', svd='numpy_svd',
                               tol=10e-7, random_state=None, verbose=0):
    """Original: tensorly.decomposition.non_negative_parafac ver 0.4.3"""
    epsilon = 10e-12

    nn_factors = initialize_factors(tensor, rank, init=init, svd=svd,
                                    random_state=random_state, non_negative=True)

    n_factors = len(nn_factors)
    norm_tensor = tl.norm(tensor, 2)
    rec_errors = []

    for iteration in range(n_iter_max):
        for mode in range(tl.ndim(tensor)):
            # khatri_rao(factors).tl.dot(khatri_rao(factors))
            # simplifies to multiplications
            sub_indices = [i for i in range(n_factors) if i != mode]
            for i, e in enumerate(sub_indices):
                if i:
                    accum = accum * tl.dot(tl.transpose(nn_factors[e]), nn_factors[e])
                else:
                    accum = tl.dot(tl.transpose(nn_factors[e]), nn_factors[e])

            numerator = tl.dot(unfold(tensor, mode), khatri_rao(nn_factors, skip_matrix=mode))
            numerator = tl.clip(numerator, a_min=epsilon, a_max=None)
            denominator = tl.dot(nn_factors[mode], accum)
            denominator = tl.clip(denominator, a_min=epsilon, a_max=None)
            nn_factors[mode] = nn_factors[mode] * numerator / denominator

        rec_error = tl.norm(tensor - kruskal_to_tensor(nn_factors), 2) / norm_tensor
        rec_errors.append(rec_error)
        if iteration > 1 and verbose:
            print('reconstruction error={}, variation={}.'.format(
                rec_errors[-1], rec_errors[-2] - rec_errors[-1]))

        if iteration > 1 and abs(rec_errors[-2] - rec_errors[-1]) < tol:
            if verbose:
                print('converged in {} iterations.'.format(iteration))
            break

    # MODIF: Keep track of final R^2 and no. of iterations
    return nn_factors, rec_error, iteration + 1


def write_decomposition_measures_to_logfile(measures, args, logfile):
    log_filename = os.path.join(args.output_dir, logfile)
    if os.path.exists(log_filename):
        write_flag = 'a'
    else:
        write_flag = 'w'

    with open(log_filename, write_flag) as metrics_file:
        if write_flag == 'w':
            metrics_file.write('ncomp,Rsq,iterations\n')

        for decomp_measure in measures:
            if decomp_measure is not None:
                n_components, rsq, n_iters = decomp_measure
                metrics_file.write(f'{n_components},{rsq},{n_iters}\n')


def get_args():
    desc = ('Wrapper script that will create PARAFAC decomposition models'
            ' for 1 to F components')
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i',
                        '--input_tensor_file',
                        required=True,
                        type=str,
                        help='NPY file holding tensor to decompose')
    parser.add_argument('-p',
                        '--tensor_properties',
                        required=True,
                        type=str,
                        help='JSON file with tensor properties not captured as numpy array')
    parser.add_argument('-f',
                        '--min_components',
                        required=False,
                        type=int,
                        default=2,
                        help='Min number of components to decompose for. '
                             '(default: %(default)s)')
    parser.add_argument('-F',
                        '--max_components',
                        required=False,
                        type=int,
                        default=20,
                        help='Max number of components to decompose for. '
                             '(default: %(default)s)')
    parser.add_argument('--init',
                        required=False,
                        type=str,
                        default='random',
                        help='PARAFAC initialization method: random, svd'
                             ' (default: %(default)s)')
    parser.add_argument('-r',
                        '--max_iter',
                        required=False,
                        type=int,
                        default=2000,
                        help='Max number of iterations when fitting PARAFAC.'
                             ' (default: %(default)s)')
    parser.add_argument('-t',
                        '--tolerance',
                        required=False,
                        type=float,
                        default=1e-7,
                        help='Tolerance when fitting PARAFAC (default: %(default)s)')
    parser.add_argument('-g',
                        '--use_gpu',
                        action='store_true',
                        help='Whether the decomposition should use GPUs')
    parser.add_argument('--pipe_id',
                        required=False,
                        type=str,
                        default='',
                        help='Unique ID of the MPS pipe to use. Default:  %(default)s)')
    parser.add_argument('--avail_ram_gb',
                        required=False,
                        type=int,
                        default=8,
                        help='Total RAM available on device (CPU or GPU)')
    parser.add_argument('--backend',
                        required=False,
                        type=str,
                        default='numpy',
                        help='Tensor backend to use')
    parser.add_argument('--skip_completed',
                        action='store_true',
                        help='Skip decomposition if mdoel files exist.')

    if len(sys.argv) > 1 and sys.argv[1].startswith('--test'):
        if sys.argv[1] == '--test_cluster':
            args = get_test_args_cluster()
        elif sys.argv[1] == '--test':
            args = get_test_args()
        logger.warning('Running on test slice')
        logger.info('TEST args:\n%s' % str(args))
    else:
        args = parser.parse_args()
        args.output_dir = os.path.dirname(args.input_tensor_file) or '.'
        args.device_id = get_device_identifier(args)
        args.rnd_seed = 123
        args.component_range = list(range(args.min_components, args.max_components + 1))

    if args.backend != 'pytorch':
        raise NotImplementedError

    return args


def get_test_args():
    args = argparse.Namespace(
        input_tensor_file='test/decomposition/slices/swath_lower_adjusted=751.50/rt_window=30.0/slice_tensor.npy',
        tensor_properties='test/decomposition/slices/swath_lower_adjusted=751.50/rt_window=30.0/tensor_properties.json',
        min_components=25,
        max_components=25,
        init='random',
        use_gpu=False,
        max_iter=5000,
        tolerance=1e-7,
        backend='pytorch',
        avail_ram_gb=16,
        device_number=0,
        skip_completed=False)
    args.output_dir = os.path.dirname(args.input_tensor_file) or '.'
    args.device_id = get_device_identifier(args)
    args.rnd_seed = 123
    args.component_range = list(range(args.min_components, args.max_components + 1))
    return args


def get_test_args_cluster():
    args = argparse.Namespace(
        input_tensor_file='test/decomposition/slices/swath_lower_adjusted=751.50/rt_window=30.0/slice_tensor.npy',
        tensor_properties='test/decomposition/slices/swath_lower_adjusted=751.50/rt_window=30.0/tensor_properties.json',
        min_components=15,
        max_components=64,
        init='random',
        use_gpu=True,
        max_iter=5000,
        tolerance=1e-7,
        backend='pytorch',
        avail_ram_gb=32,
        device_number=0,
        skip_completed=False)
    args.output_dir = os.path.dirname(args.input_tensor_file) or '.'
    args.device_id = get_device_identifier(args)
    args.rnd_seed = 123
    args.component_range = list(range(args.min_components, args.max_components + 1))
    return args


def test_decompose():
    from astropy.convolution import Gaussian1DKernel

    rnd_seed = 123
    np.random.seed(rnd_seed)

    spectra = np.array([
        [ 0,  0],
        [ 2,  4],
        [ 2,  2],
        [ 0, 10],
        [10,  0]
    ])

    expected_mass_mode = np.array([
        [0., 0.],
        [0.96, 0.49],
        [0.76, 0.76],
        [1.25, 0.0098],
        [0.017, 1.32]], dtype=np.float32)

    RT_LENGTH = 20
    NSAMPLES = 100

    elution_profile = Gaussian1DKernel(stddev=1, x_size=RT_LENGTH)
    elution_profile = elution_profile.array.reshape((1, RT_LENGTH))

    tensor = []
    for i in range(NSAMPLES):
        spectrum_mix = np.random.rand(spectra.shape[1], 1)
        scale = np.random.rand()
        sample_mixing_matrix = np.dot(spectrum_mix, elution_profile * scale)

        sample_map = np.dot(spectra, sample_mixing_matrix).T
        tensor.append(sample_map)
    tensor = np.stack(tensor)

    tensor, _ = preprocess(tensor)

    model, reconcstr_error, _ = modif_non_negative_parafac(tensor,
                                                           rank=spectra.shape[1],
                                                           n_iter_max=5000,
                                                           init='random',
                                                           tol=1e-7,
                                                           random_state=rnd_seed,
                                                           verbose=False)
    modes = [np.array(mode) for mode in model]

    assert modes[0].shape == (100, 2)
    assert modes[1].shape == (20, 2)
    assert modes[2].shape == (5, 2)

    mass_mode = np.copy(modes[2])
    mass_mode[mass_mode <= 1e-5] = 0
    # Equal to 2 decimal places is good enough
    # (randomized results might differ between machines regardless of seed)
    assert np.allclose(mass_mode, expected_mass_mode, atol=1e-2)


if __name__ == '__main__':
    main()
