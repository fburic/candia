import argparse
import csv
import inspect
import os
import sys

# _this_filename = inspect.getframeinfo(inspect.currentframe()).filename
# _this_path = os.path.dirname(os.path.abspath(_this_filename))
# sys.path.append(os.path.join(_this_path, '..'))

import msproc


def main():
    args = get_args()
    if not os.path.exists(os.path.dirname(args.output_file)):
        os.makedirs(os.path.dirname(args.output_file))

    if args.swath_isolation_windows_file is not None:
        with open(args.swath_isolation_windows_file, 'r') as windows_file:
            isolation_windows = []
            for row in csv.reader(windows_file):
                try:
                    isolation_windows.append((float(row[0]), float(row[1])))
                except ValueError:
                    # Header present
                    pass

    msproc.mzml_to_csv(mzml_filename=args.input_file,
                       swath_isolation_windows=isolation_windows,
                       min_intensity=args.min_intensity,
                       csv_filename=args.output_file,
                       overwrite=True,
                       csv_buffer_size=args.buffer_size)


def create_isolation_windows_file(mzml_filename, filename):
    isol_windows = msproc.get_swath_intervals_from_mzml(mzml_filename)
    msproc.save_isolation_windows_to_csv(isol_windows, filename)


def get_args():
    parser = argparse.ArgumentParser(description="Convert mzML to CSV")
    parser.add_argument('-i',
                        '--input_file',
                        required=True,
                        type=str,
                        help='Input samples directory')
    parser.add_argument('-s',
                        '--swath_isolation_windows_file',
                        required=False,
                        type=str,
                        help='CSV with SWATH isolation m/z windows')
    parser.add_argument('-f',
                        '--min_intensity',
                        required=False,
                        type=float,
                        default=0.0,
                        help='Filter out points with intensity lower than this value'
                             ' (default: %(default)s)')
    parser.add_argument('-o',
                        '--output_file',
                        required=True,
                        type=str,
                        help='output file ')
    parser.add_argument('-b',
                        '--buffer_size',
                        required=False,
                        type=int,
                        default=5000,
                        help='Number of CSV entries to buffer in memory (default: %(default)s)')
    return parser.parse_args()


if __name__ == '__main__':
    main()
