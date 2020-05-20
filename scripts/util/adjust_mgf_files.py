"""
Add headers required by DeepNovo

DeepNovo parsing:
https://github.com/nh2tran/DeepNovo/blob/9e31d54d31deb66ff6c13de3495754f462624e15/deepnovo_worker_io.py#L214

MGF specification:
http://www.matrixscience.com/help/data_file_help.html
"""
import argparse
import copy
import inspect
import logging
import sys
from pathlib import Path

import coloredlogs
import pyteomics.mgf as mgf
from tqdm import tqdm

_this_filename = inspect.getframeinfo(inspect.currentframe()).filename
_this_path = Path(_this_filename).parent.resolve()
sys.path.append(str(_this_path.parent))

import msproc

logging.basicConfig(format=msproc.LOG_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)
coloredlogs.install(fmt=msproc.LOG_FORMAT, level='INFO', logger=logger)


def main():
    args = get_args()

    header = mgf.read_header(args.input)
    adjusted_spectra = []

    logger.info('Reading input MGF ...')
    with mgf.read(args.input, use_header=False, read_charges=False) as reader:
        i = 0
        for spectrum in tqdm(reader):
            adjust_spectrum_metadata(args, spectrum, current_spectrum_num=i, charge=2)
            adjusted_spectra.append(spectrum)
            i += 1

    with mgf.read(args.input, use_header=False, read_charges=False) as reader:
        i = 0
        for spectrum in tqdm(reader):
            adjust_spectrum_metadata(args, spectrum, current_spectrum_num=i, charge=3)
            adjusted_spectra.append(spectrum)
            i += 1

    logger.info('Writing output MGF ...')
    mgf.write(spectra=adjusted_spectra, header='', output=args.output)
    logger.info('done.')


def adjust_spectrum_metadata(args, spectrum, current_spectrum_num, charge=2):
    spectrum_num = spectrum['params']['title'].split('scan=')[1].split('_')[0]
    spectrum['params']['title'] = spectrum_num
    spectrum['params']['scans'] = spectrum_num
    spectrum['params']['charge'] = charge
    if args.sequencer == 'deepnovo':
        spectrum['params']['rtinseconds'] = 10 + current_spectrum_num
        spectrum['params']['seq'] = 'PEPTIDE'
    return spectrum


def get_args():
    desc = 'Adjust MGF file entries'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i',
                        '--input',
                        required=True,
                        type=str,
                        help='MGF scan file from PARAFAC')
    parser.add_argument('-o',
                        '--output',
                        required=True,
                        type=str,
                        help='Adjusted MGF scan file')
    parser.add_argument('--sequencer',
                        required=True,
                        type=str,
                        help='Options: novor, deepnovo')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()