import argparse
import json
import logging
import os
from glob import glob

import coloredlogs
import numpy as np
import pandas as pd

LOG_FORMAT = "[%(asctime)s] [PID %(process)d] %(levelname)s:\t%(filename)s:%(funcName)s():%(lineno)s:\t%(message)s"
logging.basicConfig(format=LOG_FORMAT)
logger = logging.getLogger(__name__)
coloredlogs.install(fmt=LOG_FORMAT, level='INFO', logger=logger)


def main():
    args = get_args()
    slice_maps, sample_filenames = read_slice_map(args)

    try:
        slice_tensor, tensor_properties = tensorize_slice(args, slice_maps)

    except (IndexError, EmptyMSLevelException, TooFewPointsAcrossSamplesException) as e:
        logger.warning('Error while tensorizing: %s' % str(e))
        slice_tensor = np.zeros((1, 1, 1))
        tensor_properties = {'error': 'Could not tensorize'}

    write_tensor(slice_tensor, args.slice_dir)
    tensor_properties['samples'] = sample_filenames
    save_tensor_properties(tensor_properties, args.slice_dir)


def tensorize_slice(args, slice_maps):
    """
    For each sample CSV:
        - filter
        - bin RT values in cycles
        - aggregate m/z according to scan resolution

    Concatenate samples, then pivot to create data structure as:
    (sample x scan cycle x m/z partition)

    Expected columns in slice CSVs:
    ['spectrum_index', 'level', 'rt', 'mz', 'intensity', 'prec_mz', 'swath_upper_adjusted']
    """
    logger.info("Tensorizing slice %s" % args.slice_dir)

    slice_maps = slice_maps.groupby('sample_no', sort=False).apply(bin_rt_into_scan_cycles)
    slice_maps = slice_maps.drop(columns=['spectrum_index', 'rt', 'swath_upper_adjusted'])
    slice_maps = partition_mz_values(slice_maps, mz_tol=args.mz_tolerance)

    # Keep track of MS level in m/z indices + turn partition indices into str for agg
    slice_maps['mz_partition_start'] = (slice_maps['mz_partition_start']
                                        + '_ms' + slice_maps['level'].astype(str))
    slice_maps = slice_maps.drop(columns=['level'])

    slice_maps = aggregate_mz_partitions_in_concatenated_sample_maps(slice_maps,
                                                                     partition_agg_func=np.nansum)
    slice_maps = remove_infrequent_points(slice_maps)
    if slice_maps.shape[0] == 0:
        raise TooFewPointsAcrossSamplesException(
            'Too few data points across all samples. Skipping slice.'
        )

    # Align m/z indices by pivoting on concatenated maps (include all m/z values)
    slice_maps = slice_maps.pivot_table(index=['sample_no', 'cycle'],
                                        columns='mz_partition_start',
                                        values='intensity',
                                        aggfunc=np.nansum)

    # Sort m/z indices (columns) by level and value
    slice_maps = slice_maps.reindex(sorted(slice_maps.columns,
                                           key=lambda x: (int(x.split('_ms')[1]),
                                                          float(x.split('_ms')[0]))),
                                    axis=1)
    mz_indices = slice_maps.columns.values.tolist()
    slice_maps = list(slice_maps.groupby('sample_no').apply(lambda df: df.values))

    # Needed for padding the resulting matrices to same size
    max_rt_len = max([sample_map.shape[0] for sample_map in slice_maps])

    slice_tensor = np.zeros(shape=(len(slice_maps), max_rt_len, len(mz_indices)),
                            dtype=np.float32)
    slice_tensor.fill(np.nan)

    # Convert pandas.DataFrames to numpy.arrays
    for i in range(len(slice_maps)):
        len_time = slice_maps[i].shape[0]
        len_mz = slice_maps[i].shape[1]

        slice_tensor[i][0:len_time, 0:len_mz] = slice_maps[i]

    tensor_properties = {'mz_indices': mz_indices}

    return slice_tensor, tensor_properties


def bin_rt_into_scan_cycles(sample_map):
    """
    Bin the rt values into cycles so that
    each MS 1 scan is on the same time axis as its associated MS 2 scans

    The binning is done using the acquisition times of MS 1 scans so that spectra
    are grouped as (MS 1, MS 2, MS 2, ..), (MS 1, MS 2, MS 2, ...)

    :return:
    sample_map =
    The sample map with a new colum "cycle" marking each entry with
    the scan cycle it belongs to

    ms1_acquisition_times =
    Acquisition times of MS 1 the bins used to bin spectra into cycles,
    plus an extra trailing time used to close the last bin.
    """
    ms1_acquisition_times = sorted(
        sample_map[sample_map['level'] == 1]['rt'].unique()
    )
    if len(ms1_acquisition_times) == 0:
        raise EmptyMSLevelException('No MS 1 acquisition times found')

    # Complete last bin, but make sure to include last value
    # since since working with [closed, open) intervals
    ms1_acquisition_times.append(ms1_acquisition_times[-1] + 0.1)

    map_binned = sample_map.assign(
        cycle = pd.cut(sample_map['rt'],
                       bins=ms1_acquisition_times,
                       right=False,
                       include_lowest=True,
                       labels=range(len(ms1_acquisition_times) - 1)))

    # These are MS 2 spectra ("tails") that were cut from the neighbouring slice
    num_nans = len(map_binned[map_binned['cycle'].isna()]['spectrum_index'].unique())
    if num_nans > 0:
        n_spectra = len(map_binned['spectrum_index'].unique())
        logger.warning(
            'Sample %s:\tDropping %d MS 2 spectra (frac = %f) from neighbouring slice' %
            (map_binned['file'].values[0], num_nans, 100 * num_nans / n_spectra))

    # Drop MS 2 spectra ("tails") that were cut from the neighbouring slice
    map_binned = map_binned.dropna()
    map_binned['cycle'] = map_binned['cycle'].astype(int)  # uniform treatment

    return map_binned


def filter(sample_map, intensity_cutoff):
    sample_map = sample_map[sample_map['intensity'] >= intensity_cutoff]
    logger.info('Filtered out intensities < ' + str(intensity_cutoff))
    return sample_map


def remove_infrequent_points(slice_maps):
    min_n_points_per_sample = 5

    n_points = slice_maps.groupby(['sample_no', 'mz_partition_start'])['cycle'].count()
    most_points_across_samples = n_points.groupby('mz_partition_start').max()
    mz_with_enough_points = most_points_across_samples[
        most_points_across_samples >= min_n_points_per_sample
        ]
    mz_with_enough_points = set(mz_with_enough_points.index.values.tolist())

    filtererd_slice_maps = slice_maps[
        slice_maps['mz_partition_start'].isin(mz_with_enough_points)
    ]
    return filtererd_slice_maps


def aggregate_mz_partitions_in_concatenated_sample_maps(sample_maps,
                                                        partition_agg_func=np.nansum):
    # Temporarily replace any NaNs since groupby() drops them
    sample_maps.fillna(-1, inplace=True)
    sample_map = sample_maps.groupby(['sample_no', 'cycle', 'mz_partition_start'],
                                      as_index=False).agg(partition_agg_func)

    sample_map.replace(-1, np.nan, inplace=True)
    return sample_map


def partition_mz_values(msmap, mz_tol):
    """
    Return msmap with a new column 'mz_partition_start' which holds
    the min m/z value of the partition to which each m/z value belongs, as str
    e.g.
    mz = 100.01  mz_partition_start = "100.01"
    mz = 100.02  mz_partition_start = "100.01"
    mz = 100.20  mz_partition_start = "100.20"
    """
    partitioned_map = []
    for level in [1, 2]:
        msmap_level = msmap[msmap['level'] == level]
        if msmap_level.shape[0] == 0:
            raise EmptyMSLevelException(f'MS {level} completely empty')

        mz_partitions = get_mz_partitions(sorted(msmap_level['mz'].unique()),
                                          mz_tol)
        msmap_level = msmap_level.assign(
            mz_partition_start=
                pd.cut(msmap_level['mz'],
                       bins=mz_partitions,
                       right=False,
                       include_lowest=True,
                       labels=mz_partitions[:-1])
                .apply(lambda mz: "{:.4f}".format(mz))
        )
        partitioned_map.append(msmap_level)

    partitioned_map = pd.concat(partitioned_map, ignore_index=True)
    return partitioned_map


def get_mz_partitions(mz_list, mz_tol):
    # Create list of m/z bins.
    # Each bin (i.e. m/z partition) will have its width = ppm(first m/z)
    mz_partitions = [min(mz_list)]
    current_mz_partition_start = mz_partitions[-1]
    len_partition = ppm_tol(current_mz_partition_start, mz_tol)

    for i, mz in enumerate(mz_list):
        dist_from_partition_start = mz - current_mz_partition_start
        if dist_from_partition_start > len_partition:
            mz_partitions.append(mz)
            current_mz_partition_start = mz_partitions[-1]
            len_partition = ppm_tol(current_mz_partition_start, mz_tol)

    # Close interval to be able to include last value
    mz_partitions.append(mz_partitions[-1] + 1)
    return mz_partitions


def ppm_tol(mz, ppm=40):
    return (mz * ppm) / float(1e6)


def get_sample_directories(args):
    sample_filename_pattern = 'file=*'
    sample_directories = glob(os.path.join(args.slice_dir, sample_filename_pattern))
    if len(sample_directories) == 0:
        msg = "Slice %s has no samples. Aborting" % args.slice_dir
        logger.error(msg)
        raise Exception(msg)
    return sample_directories


def read_slice_map(args):
    slice_maps = pd.read_csv(glob(os.path.join(args.slice_dir, 'part-*.csv'))[0])

    sample_filenames = sorted(list(slice_maps['file'].unique()))
    sample_file_to_num = {filename: i for i, filename in enumerate(sample_filenames)}
    slice_maps = slice_maps.assign(
        sample_no = slice_maps['file'].apply(lambda fname: sample_file_to_num[fname]
    ))
    return slice_maps, sample_filenames


def read_csv_parts(sample_dir):
    csv_partitions = glob(os.path.join(sample_dir, "*.csv"))
    if not csv_partitions:
        raise Exception("Sample %s has no data. Aborting" % sample_dir)

    sample_map = pd.concat((pd.read_csv(f) for f in csv_partitions),
                           ignore_index=True)
    return sample_map


def write_tensor(slice_tensor, slice_dir):
    tensor_path = os.path.join(slice_dir, 'slice_tensor.npy')
    np.save(tensor_path, slice_tensor)
    logger.info('Tensor shape: %s' % str(slice_tensor.shape))
    logger.info('Wrote %s' % tensor_path)


def save_tensor_properties(tensor_properties, slice_dir):
    with open(os.path.join(slice_dir, 'tensor_properties.json'), 'w') as properties_file:
        json.dump(tensor_properties, properties_file, indent=4)


class EmptyMSLevelException(Exception):
    pass


class TooFewPointsAcrossSamplesException(Exception):
    pass


def get_args():
    desc = ("Convert (swath, rt) slice generated by Spark splitting script "
            "to tensors as .nyp files.")
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-s',
                        '--slice_dir',
                        required=True,
                        type=str,
                        help='Single slice directory organized as swath/rt_window. '
                             'Output tensor will be written here')
    parser.add_argument('-d',
                        '--mz_decimals',
                        required=False,
                        type=int,
                        default=10,
                        help="Round down m/z values to this many decimals (default: %(default)s) "
                             "while summing intensity values")
    parser.add_argument('-c',
                        '--intensity_cutoff',
                        required=False,
                        type=float,
                        help="Keep only intensity values higher or equal than this")

    parser.add_argument('-t',
                        '--mz_tolerance',
                        required=False,
                        type=float,
                        default=40,
                        help="m/z tolerance in ppm (default: %(default)s)")

    return parser.parse_args()


# FIXME: Update to single input CSV / slice
# def test_tensorize_slice():
#     column_labels = ['spectrum_index', 'rt', 'mz', 'intensity', 'level', 'swath_upper_adjusted']
#
#     sample_11 = pd.DataFrame([
#         [123, 1.01, 10, 100, 1, 0],
#         [124, 1.02, 20, 200, 2, 0]
#     ], columns=column_labels)
#
#     sample_12 = pd.DataFrame([
#         [345, 3.02, 30, 300, 1, 0],
#         [346, 3.03, 30, 333, 2, 0]
#     ], columns=column_labels)
#
#     sample_2 = pd.DataFrame([
#         [347, 2.01, 11, 400, 1, 0],
#         [348, 2.06, 22, 500, 2, 0],
#         [349, 2.07, 22, 800, 1, 0],
#         [350, 2.08, 33, 600, 2, 0]
#     ], columns=column_labels)
#
#     expected_mz_indices = ['10.0000000000_ms1',
#                            '11.0000000000_ms1',
#                            '22.0000000000_ms1',
#                            '30.0000000000_ms1',
#                            '20.0000000000_ms2',
#                            '22.0000000000_ms2',
#                            '30.0000000000_ms2',
#                            '33.0000000000_ms2']
#
#     expected_tensor = np.array([
#         [[np.nan,   400., np.nan, np.nan, np.nan,   500., np.nan, np.nan],
#          [np.nan, np.nan,   800., np.nan, np.nan, np.nan, np.nan,   600.]],
#         [[  100., np.nan, np.nan, np.nan,   200., np.nan, np.nan, np.nan],
#          [np.nan, np.nan, np.nan,   300., np.nan, np.nan,   333., np.nan]]],
#     dtype=np.float32)
#
#     expected_tensor = expected_tensor.astype(np.float32)
#
#     test_dir = '/tmp/test_tensorize'
#     if not os.path.exists(test_dir + '/file=1'):
#         os.makedirs(test_dir + '/file=1')
#     if not os.path.exists(test_dir + '/file=2'):
#         os.makedirs(test_dir + '/file=2')
#
#     sample_11.to_csv('/tmp/test_tensorize/file=1/part-1.csv', index=False)
#     sample_12.to_csv('/tmp/test_tensorize/file=1/part-2.csv', index=False)
#     sample_2.to_csv('/tmp/test_tensorize/file=2/part-1.csv', index=False)
#
#     args = argparse.Namespace(
#         slice_dir=test_dir,
#         mz_decimals=10,
#         intensity_cutoff=0,
#         mz_tolerance=40
#     )
#
#     slice_maps, sample_filenames = read_and_concatenate_slice_sample_maps(args)
#     slice_tensor, tensor_properties = tensorize_slice(args, slice_maps)
#
#     assert tensor_properties['mz_indices'] == expected_mz_indices
#     assert np.allclose(slice_tensor, expected_tensor, equal_nan=True)


def test_mz_partitioning():
    from io import StringIO

    test_map_csv = StringIO("""
level,rt,mz,swath_upper_adjusted,sample_no,spectrum_index,intensity
2,2646.76806640625,350.1284179688,474.5,39,26974,1.31461501121521
2,2667.387939453125,350.12890625,474.5,70,27187,1.509950876235962
2,2643.34912109375,350.12890625,474.5,73,26940,1.5772640705108645
2,2657.08203125,350.1294555664,474.5,80,27078,1.643269658088684
2,2643.3359375,350.1296081543,474.5,39,26939,2.493281841278076
2,2663.943115234375,350.129699707,474.5,60,27150,1.5086196660995483
2,2657.10400390625,350.1301269531,474.5,87,27079,1.970595121383667
2,2650.2109375,350.1304626465,474.5,9,27003,1.5521864891052246
2,2667.373046875,350.1306152344,474.5,55,27172,3.075400352478028
2,2670.821044921875,350.1307983398,474.5,70,27222,1.776573657989502
2,2650.200927734375,350.1309814453,474.5,74,27012,3.3829162120819087
2,2653.64306640625,350.1317749023,474.5,9,27038,1.1958253383636477
2,2653.634033203125,350.1321411133,474.5,54,27042,2.0933663845062256
2,2650.201904296875,350.1328735352,474.5,78,26970,1.3779486417770386
2,2643.34912109375,350.1332702637,474.5,82,26939,1.4202085733413696
2,2643.35400390625,350.1333312988,474.5,34,26939,4.3556809425353995
2,2646.777099609375,350.1335449219,474.5,9,26968,1.0941891670227049
2,2663.94189453125,350.1338195801,474.5,84,27151,2.3810997009277344
2,2646.7900390625,350.1339111328,474.5,69,26977,1.871619462966919
2,2657.06396484375,350.1363220215,474.5,74,27082,1.0123751163482666
2,2667.389892578125,350.1369628906,474.5,76,27180,1.1245584487915041
2,2646.76904296875,350.1375732422,474.5,78,26935,2.0114877223968506
2,2650.215087890625,350.1397399902,474.5,30,27006,1.3119505643844604
2,2653.6708984375,350.1416625977,474.5,87,27044,1.1535953283309937
2,2663.946044921875,350.1422119141,474.5,80,27148,1.0546185970306396
2,2646.781982421875,350.1463012695,474.5,84,26976,1.0546137094497678
2,2653.633056640625,350.1466064453,474.5,42,27012,1.1668118238449097
2,2643.35888671875,350.1487731934,474.5,69,26942,1.2526088953018188
2,2653.657958984375,350.1489257812,474.5,70,27047,1.3199162483215332
2,2698.251953125,350.1536254883,474.5,43,27467,1.0546122789382937
2,2670.806884765625,350.1537780762,474.5,84,27221,1.2011445760726929
2,2646.77294921875,350.154510498,474.5,86,26968,1.0717763900756836
2,2684.531982421875,350.1552124023,474.5,55,27347,1.4360883235931396
2,2650.2080078125,350.1553344727,474.5,35,27009,1.2275344133377075
2,2674.201904296875,350.1555786133,474.5,37,27251,1.0242842435836792
2,2650.218994140625,350.1562194824,474.5,28,26982,1.1430357694625854
2,2674.22900390625,350.1563110352,474.5,78,27215,1.0334919691085815
2,2674.24609375,350.1563415527,474.5,38,27259,1.847889423370361
2,2663.94091796875,350.1564025879,474.5,55,27137,1.2697982788085938
2,2646.779052734375,350.1564025879,474.5,68,26945,1.0150121450424194
2,2687.972900390625,1249.6260986328,474.5,28,27367,2.1194310188293457
2,2698.26904296875,1249.6295166016,474.5,56,27451,1.835181474685669
2,2681.094970703125,1249.6326904297,474.5,19,27287,1.1769421100616455
2,2687.969970703125,1249.6334228516,474.5,56,27346,1.2966029644012451
2,2691.406005859375,1249.6339111328,474.5,28,27402,1.0397763252258299
2,2667.37109375,1249.6346435547,474.5,81,27144,1.8252184391021729
2,2657.0830078125,1249.6357421875,474.5,49,27030,1.7828432321548462
2,2657.073974609375,1249.6363525391,474.5,83,27033,1.2841426134109497
2,2684.529052734375,1249.6384277344,474.5,31,27306,1.588356852531433
2,2694.837890625,1249.638671875,474.5,48,27436,1.2317742109298706
2,2687.968017578125,1249.6416015625,474.5,51,27364,1.9199864864349363
2,2687.9619140625,1249.6446533203,474.5,12,27359,2.0097687244415283
2,2646.778076171875,1249.6450195312,474.5,90,26948,1.615800857543945
2,2670.8310546875,1249.6458740234,474.5,27,27218,1.0223441123962402
2,2687.967041015625,1249.6469726562,474.5,59,27359,1.80779230594635
2,2684.529052734375,1249.6506347656,474.5,12,27324,2.884999513626098
2,2646.781982421875,1249.6511230469,474.5,89,26940,1.7479568719863892
2,2677.657958984375,1249.6512451172,474.5,42,27257,1.6756658554077148
2,2691.40087890625,1249.6516113281,474.5,51,27399,1.2642005681991575
2,2687.972900390625,1249.6574707031,474.5,48,27366,1.5933473110198977
2,2670.798095703125,1249.6616210938,474.5,12,27184,1.2741961479187012
2,2681.10400390625,1249.6619873047,474.5,88,27299,1.755459189414978
2,2674.24609375,1249.6639404297,474.5,64,27228,1.037308931350708
2,2663.950927734375,1249.6680908203,474.5,62,27124,1.9325196743011477
2,2694.8359375,1249.6763916016,474.5,23,27433,1.19189453125
2,2698.251953125,1249.6834716797,474.5,39,27499,1.0123877525329592
2,2677.676025390625,1249.6843261719,474.5,28,27262,1.0622352361679075
2,2687.968017578125,1249.6953125,474.5,51,27364,1.039810061454773
2,2687.972900390625,1250.6291503906,474.5,48,27366,1.3270496129989624
2,2674.239013671875,1250.6375732422,474.5,59,27219,1.4817343950271606
2,2646.782958984375,1250.6427001953,474.5,77,26947,1.6189217567443848
2,2670.803955078125,1250.6428222656,474.5,81,27179,1.7236827611923218
2,2663.943115234375,1250.6431884766,474.5,61,27132,1.1025760173797607
2,2667.37109375,1250.6500244141,474.5,81,27144,1.1025571823120115
2,2684.531982421875,1250.6549072266,474.5,35,27359,1.0177780389785769
2,2694.8359375,1250.6943359375,474.5,11,27461,1.2946630716323853
2,2687.9619140625,1250.6953125,474.5,12,27359,1.2397948503494265
2,2667.3759765625,1250.6976318359,474.5,61,27167,1.1225671768188477
2,2667.375,1267.6629638672,474.5,88,27159,1.1402027606964111
2,2677.696044921875,1268.7048339844,474.5,27,27288,1.0351276397705078
""")

    expected_map_csv = StringIO("""
mz_partition_start,level,rt,mz,swath_upper_adjusted,sample_no,spectrum_index,intensity
350.1284179688,2,2646.76806640625,350.1284179688,474.5,39,26974,1.31461501121521
350.1284179688,2,2667.387939453125,350.12890625,474.5,70,27187,1.509950876235962
350.1284179688,2,2643.34912109375,350.12890625,474.5,73,26940,1.5772640705108645
350.1284179688,2,2657.08203125,350.1294555664,474.5,80,27078,1.643269658088684
350.1284179688,2,2643.3359375,350.1296081543,474.5,39,26939,2.493281841278076
350.1284179688,2,2663.943115234375,350.129699707,474.5,60,27150,1.5086196660995483
350.1284179688,2,2657.10400390625,350.1301269531,474.5,87,27079,1.970595121383667
350.1284179688,2,2650.2109375,350.1304626465,474.5,9,27003,1.5521864891052246
350.1284179688,2,2667.373046875,350.1306152344,474.5,55,27172,3.075400352478028
350.1284179688,2,2670.821044921875,350.1307983398,474.5,70,27222,1.776573657989502
350.1284179688,2,2650.200927734375,350.1309814453,474.5,74,27012,3.3829162120819087
350.1284179688,2,2653.64306640625,350.1317749023,474.5,9,27038,1.1958253383636477
350.1284179688,2,2653.634033203125,350.1321411133,474.5,54,27042,2.0933663845062256
350.1284179688,2,2650.201904296875,350.1328735352,474.5,78,26970,1.3779486417770386
350.1284179688,2,2643.34912109375,350.1332702637,474.5,82,26939,1.4202085733413696
350.1284179688,2,2643.35400390625,350.1333312988,474.5,34,26939,4.3556809425353995
350.1284179688,2,2646.777099609375,350.1335449219,474.5,9,26968,1.0941891670227049
350.1284179688,2,2663.94189453125,350.1338195801,474.5,84,27151,2.3810997009277344
350.1284179688,2,2646.7900390625,350.1339111328,474.5,69,26977,1.871619462966919
350.1284179688,2,2657.06396484375,350.1363220215,474.5,74,27082,1.0123751163482666
350.1284179688,2,2667.389892578125,350.1369628906,474.5,76,27180,1.1245584487915041
350.1284179688,2,2646.76904296875,350.1375732422,474.5,78,26935,2.0114877223968506
350.1284179688,2,2650.215087890625,350.1397399902,474.5,30,27006,1.3119505643844604
350.1284179688,2,2653.6708984375,350.1416625977,474.5,87,27044,1.1535953283309937
350.1284179688,2,2663.946044921875,350.1422119141,474.5,80,27148,1.0546185970306396
350.1463012695,2,2646.781982421875,350.1463012695,474.5,84,26976,1.0546137094497678
350.1463012695,2,2653.633056640625,350.1466064453,474.5,42,27012,1.1668118238449097
350.1463012695,2,2643.35888671875,350.1487731934,474.5,69,26942,1.2526088953018188
350.1463012695,2,2653.657958984375,350.1489257812,474.5,70,27047,1.3199162483215332
350.1463012695,2,2698.251953125,350.1536254883,474.5,43,27467,1.0546122789382937
350.1463012695,2,2670.806884765625,350.1537780762,474.5,84,27221,1.2011445760726929
350.1463012695,2,2646.77294921875,350.154510498,474.5,86,26968,1.0717763900756836
350.1463012695,2,2684.531982421875,350.1552124023,474.5,55,27347,1.4360883235931396
350.1463012695,2,2650.2080078125,350.1553344727,474.5,35,27009,1.2275344133377075
350.1463012695,2,2674.201904296875,350.1555786133,474.5,37,27251,1.0242842435836792
350.1463012695,2,2650.218994140625,350.1562194824,474.5,28,26982,1.1430357694625854
350.1463012695,2,2674.22900390625,350.1563110352,474.5,78,27215,1.0334919691085815
350.1463012695,2,2674.24609375,350.1563415527,474.5,38,27259,1.847889423370361
350.1463012695,2,2663.94091796875,350.1564025879,474.5,55,27137,1.2697982788085938
350.1463012695,2,2646.779052734375,350.1564025879,474.5,68,26945,1.0150121450424194
1249.6247558594,2,2687.972900390625,1249.6260986328,474.5,28,27367,2.1194310188293457
1249.6247558594,2,2698.26904296875,1249.6295166016,474.5,56,27451,1.835181474685669
1249.6247558594,2,2681.094970703125,1249.6326904297,474.5,19,27287,1.1769421100616455
1249.6247558594,2,2687.969970703125,1249.6334228516,474.5,56,27346,1.2966029644012451
1249.6247558594,2,2691.406005859375,1249.6339111328,474.5,28,27402,1.0397763252258299
1249.6247558594,2,2667.37109375,1249.6346435547,474.5,81,27144,1.8252184391021729
1249.6247558594,2,2657.0830078125,1249.6357421875,474.5,49,27030,1.7828432321548462
1249.6247558594,2,2657.073974609375,1249.6363525391,474.5,83,27033,1.2841426134109497
1249.6247558594,2,2684.529052734375,1249.6384277344,474.5,31,27306,1.588356852531433
1249.6247558594,2,2694.837890625,1249.638671875,474.5,48,27436,1.2317742109298706
1249.6247558594,2,2687.968017578125,1249.6416015625,474.5,51,27364,1.9199864864349363
1249.6247558594,2,2687.9619140625,1249.6446533203,474.5,12,27359,2.0097687244415283
1249.6247558594,2,2646.778076171875,1249.6450195312,474.5,90,26948,1.615800857543945
1249.6247558594,2,2670.8310546875,1249.6458740234,474.5,27,27218,1.0223441123962402
1249.6247558594,2,2687.967041015625,1249.6469726562,474.5,59,27359,1.80779230594635
1249.6247558594,2,2684.529052734375,1249.6506347656,474.5,12,27324,2.884999513626098
1249.6247558594,2,2646.781982421875,1249.6511230469,474.5,89,26940,1.7479568719863892
1249.6247558594,2,2677.657958984375,1249.6512451172,474.5,42,27257,1.6756658554077148
1249.6247558594,2,2691.40087890625,1249.6516113281,474.5,51,27399,1.2642005681991575
1249.6247558594,2,2687.972900390625,1249.6574707031,474.5,48,27366,1.5933473110198977
1249.6247558594,2,2670.798095703125,1249.6616210938,474.5,12,27184,1.2741961479187012
1249.6247558594,2,2681.10400390625,1249.6619873047,474.5,88,27299,1.755459189414978
1249.6247558594,2,2674.24609375,1249.6639404297,474.5,64,27228,1.037308931350708
1249.6247558594,2,2663.950927734375,1249.6680908203,474.5,62,27124,1.9325196743011477
1249.6763916016,2,2694.8359375,1249.6763916016,474.5,23,27433,1.19189453125
1249.6763916016,2,2698.251953125,1249.6834716797,474.5,39,27499,1.0123877525329592
1249.6763916016,2,2677.676025390625,1249.6843261719,474.5,28,27262,1.0622352361679075
1249.6763916016,2,2687.968017578125,1249.6953125,474.5,51,27364,1.039810061454773
1250.6291503906,2,2687.972900390625,1250.6291503906,474.5,48,27366,1.3270496129989624
1250.6291503906,2,2674.239013671875,1250.6375732422,474.5,59,27219,1.4817343950271606
1250.6291503906,2,2646.782958984375,1250.6427001953,474.5,77,26947,1.6189217567443848
1250.6291503906,2,2670.803955078125,1250.6428222656,474.5,81,27179,1.7236827611923218
1250.6291503906,2,2663.943115234375,1250.6431884766,474.5,61,27132,1.1025760173797607
1250.6291503906,2,2667.37109375,1250.6500244141,474.5,81,27144,1.1025571823120115
1250.6291503906,2,2684.531982421875,1250.6549072266,474.5,35,27359,1.0177780389785769
1250.6943359375,2,2694.8359375,1250.6943359375,474.5,11,27461,1.2946630716323853
1250.6943359375,2,2687.9619140625,1250.6953125,474.5,12,27359,1.2397948503494265
1250.6943359375,2,2667.3759765625,1250.6976318359,474.5,61,27167,1.1225671768188477
1267.6629638672,2,2667.375,1267.6629638672,474.5,88,27159,1.1402027606964111
1268.7048339844,2,2677.696044921875,1268.7048339844,474.5,27,27288,1.0351276397705078
""")

    test_map = pd.read_csv(test_map_csv, float_precision='round_trip')
    expected_map = pd.read_csv(expected_map_csv, float_precision='round_trip')

    result_map = partition_mz_values(test_map, mz_tol=40)

    cols = ['mz_partition_start', 'level', 'rt', 'mz', 'swath_upper_adjusted',
            'sample_no', 'spectrum_index', 'intensity']

    np.allclose(result_map[cols].values.astype(np.float32),
                expected_map[cols].values.astype(np.float32))


if __name__ == '__main__':
    main()
