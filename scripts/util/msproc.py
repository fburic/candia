"""
Collection of mass spectrometry data processing methods.
"""
import base64
import csv
import logging
import re
from multiprocessing import Pool, cpu_count
import os
import sys

import coloredlogs
import numpy as np
import pandas as pd
import intervaltree as itree

import pymzml
from pyteomics import mzxml as pyteomics_mzxml
from pyteomics import mzml
from pyteomics import mgf
from pyteomics import mass
from tqdm import tqdm


LOG_FORMAT = "[%(asctime)s] [PID %(process)d] %(levelname)s:\t%(filename)s:%(funcName)s():%(lineno)s:\t%(message)s"
logging.basicConfig(format=LOG_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)
coloredlogs.install(fmt=LOG_FORMAT, level='INFO', logger=logger)


RT_DIGITS = 4
MZ_DIGITS = 10
I_DIGITS = 8

csv_header = ["spectrum_index",
              "level",
              "rt",
              "mz",
              "intensity",
              "prec_mz",
              "prec_isolation_window_start",
              "prec_isolation_window_end"]


adjusted_swaths_map_dtypes = {
    'level': 'category',
    'rt': np.float32,
    'mz': np.float32,
    'intensity': np.float32,
    'swath_upper_adjusted': np.float16,
    'swath_lower_adjusted': np.float16
}


MZXML_HEADER = ("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
                "<mzXML xmlns=\"http://sashimi.sourceforge.net/schema_revision/mzXML_3.2\"\n"
                "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
                "xsi:schemaLocation=\"http://sashimi.sourceforge.net/schema_revision/mzXML_3.2"
                " http://sashimi.sourceforge.net/schema_revision/mzXML_3.2/mzXML_idx_3.2.xsd\">\n"
                "<msRun scanCount=\"N/A\">\n")

np.random.seed(123)


def mzml_to_csv(mzml_filename,
                swath_isolation_windows=None,
                min_intensity=0.0,
                csv_filename=None,
                overwrite=True,
                csv_buffer_size=5000):
    """
    Retention time is saved to seconds, regardless of input file unit.

    For MS 1, the current m/z is recorded as "precursor m/z" and the
    SWATH interval in which this falls is recorded as "isolation window".
    This is done for the sake of uniformity in processing scan data.

    Form MS 2, precursor m/z = isolation window target m/z [MS:1000040]

    Note that precursor ion m/z [MS:1000744] is assumed to be equal
    to the isolation window target m/z [MS:1000827]. The former is the one saved.
    """
    if swath_isolation_windows is None:
        logger.info("No SWATH isolation windows supplied. "
                    "Getting these from the input file ...")
        isol_windows = get_swath_intervals_from_mzml(mzml_filename)
        logger.info("... done")
    else:
        isol_windows = itree.IntervalTree()
        for idx, window in enumerate(swath_isolation_windows):
            isol_windows.addi(window[0], window[1], idx)

    if csv_filename is None:
        csv_filename = mzml_filename + ".csv"

    if not overwrite and os.path.exists(csv_filename):
        logger.warning("Aborting since file exists: " + csv_filename)
        return

    time_unit = get_time_unit(mzml_filename)
    logger.info('Time unit: ' + str(time_unit))

    logger.info("Writing %s ...", csv_filename)

    processed_spectra = 0
    with mzml.read(mzml_filename) as msrun:
        with open(csv_filename, 'w') as csv_output:
            writer = csv.writer(csv_output)
            csv_contents = [csv_header]

            n_spectra_processed = 0
            n_spectra_problematic = 0

            for spectrum in tqdm(msrun, desc="Spectra processed"):
                n_spectra_processed += 1
                try:
                    spectrum_index = spectrum['index']

                    level = spectrum['ms level']
                    if level is not None:
                        level = int(level)
                    else:
                        level = ''
                        logger.warning("Spectrum index=%s missing MS level.",
                                       spectrum_index)
                        n_spectra_problematic += 1

                    if 'scan start time' not in spectrum['scanList']['scan'][0]:
                        logger.warning("Spectrum index=%s missing scan time. Skipping.",
                                       spectrum_index)
                        n_spectra_problematic += 1
                        continue

                    rt = spectrum['scanList']['scan'][0]['scan start time']
                    if time_unit == 'minute':
                        rt *= 60
                    rt = np.around(rt, decimals=RT_DIGITS)

                    if level == 2:
                        if spectrum['precursorList']['count'] == 0:
                            logger.warning("Precursor isolation window information missing "
                                           f"for spectrum {spectrum_index}. Skipping.")
                            continue

                        prec_info = spectrum['precursorList']['precursor'][0]['isolationWindow']
                        prec_isolation_window_target_mz = prec_info['isolation window target m/z']
                        prec_mz = np.around(prec_isolation_window_target_mz, decimals=MZ_DIGITS)

                    if 'm/z array' not in spectrum:
                        logger.warning(f'Spectrum {spectrum_index} has no m/z array. Skipping')
                        n_spectra_problematic += 1
                        continue

                    for mz, intensity in zip(spectrum['m/z array'], spectrum['intensity array']):

                        if level == 2 and intensity < min_intensity:
                            continue

                        mz = np.around(mz, decimals=MZ_DIGITS)
                        intensity = np.around(intensity, decimals=I_DIGITS)
                        if level == 1:
                            prec_mz = mz

                        window_containing_prec_mz = sorted(isol_windows[prec_mz])
                        if window_containing_prec_mz:
                            prec_isolation_window_start = window_containing_prec_mz[0].begin
                            prec_isolation_window_end = window_containing_prec_mz[0].end
                        else:
                            # Discard values outside of SWATH range
                            continue

                        entry = [spectrum_index,
                                 level,
                                 rt,
                                 mz,
                                 intensity,
                                 prec_mz,
                                 prec_isolation_window_start,
                                 prec_isolation_window_end]
                        csv_contents.append(entry)

                    if processed_spectra == csv_buffer_size:
                        writer.writerows(csv_contents)
                        csv_contents = []
                        processed_spectra = 0
                    else:
                        processed_spectra += 1

                except Exception as e:
                    logger.error("Exception for file: " + mzml_filename)
                    logger.error("spectrum index = %s", spectrum_index)
                    logger.error("spectrum object: %s", spectrum)
                    logger.error("SWATH intervals: %s", isol_windows)
                    os.rename(csv_filename, csv_filename + ".error")
                    csv_output.close()
                    raise e

            writer.writerows(csv_contents)
            logger.info("Wrote %s", csv_filename)
            logger.info(f"{n_spectra_problematic} / {n_spectra_processed} "
                        "spectra could not be processed")


def get_time_unit(mzml_filename):
    pattern = re.compile(re.compile('UO:0000031" unitName="(\w+)"'))
    with open(mzml_filename , 'r') as mzml_file:
        for line in mzml_file:
            match = pattern.search(line)
            if match is not None:
                return match.group(1)


def mass_mode_to_mzxml(mass_mode, mz_indices, mzxml_filename,
                       intensity_cutoff_bin=0, isolation_window_center=-1):
    """
    Take mass mode as Numpy array and mz_indices as strings,
    convert to a Pandas DataFrame indexed by the mz indices,
    then generate the corresponding indexed mxXML file.
    """
    logger.info('Wrapping mass mode in mzXML format as ' + mzxml_filename)
    mass_mode = as_dataframe_with_mz_index(mass_mode, mz_indices)

    with MassMode2MzXMLEncoder(mzxml_filename) as mzxml_encoder:
        mzxml_encoder.write_scans_from_mass_mode_df(mass_mode,
                                                    intensity_cutoff_bin=intensity_cutoff_bin,
                                                    isolation_window_center=isolation_window_center)


class MassMode2MzXMLEncoder(object):
    def __init__(self, mzxml_filename):
        self.mzxml_filename = mzxml_filename

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.finish()

    def start(self):
        with open(self.mzxml_filename, 'w') as mzxml:
            mzxml.write(MZXML_HEADER)

    def finish(self):
        with open(self.mzxml_filename, 'a') as mzxml:
            mzxml.write('</msRun>\n')
            mzxml.write('</mzXML>\n')
        record_scan_count(self.mzxml_filename)
        index_mzxml_file(self.mzxml_filename)

    def write_scans_from_mass_mode_df(self,
                                      mass_mode: pd.DataFrame,
                                      intensity_cutoff_bin=0,
                                      isolation_window_center=-1):
        """
        Converts a mass mode given DataFrame into <scan> structures,
        which are written to the output mzXML file.

        NOTE: The file is afterwards closed but not finalized
        (</mzxml> not written and the file is not indexed)

        :param mass_mode: A panda.DataFrame with columns holding
        mass mode components (spectra), as well as two columns 'mz' and 'level'
        """
        xml_content = ''
        scan_no = 0
        for component_idx in mass_mode.drop(columns=['mz', 'level']):
            spectrum = mass_mode[[component_idx, 'mz', 'level']]

            filtered_spectrum = filter_out_MS2_entries_with_intensity_in_lower_percent(spectrum,
                                                                                       intensity_cutoff_bin)
            if filtered_spectrum.empty:
                logger.warning('MS 2 data in component %d has no intensity values above threshold' % component_idx)
                xml_content += ''
            else:
                xml_spectrum_ms2 = convert_mass_mode_component_to_mzxml(filtered_spectrum,
                                                                        scan_id=component_idx,
                                                                        isolation_window_center=isolation_window_center)
                xml_content += xml_spectrum_ms2

            scan_no += 1

        with open(self.mzxml_filename, 'a') as mzxml:
            mzxml.write(xml_content)


def convert_mass_mode_component_to_mzxml(spectrum, scan_id, isolation_window_center):
    """
    Produce mzXML <scan> structure encoding an MS 2 spectrum from a single mass mode component.

    :param spectrum: A panda.DataFrame with columns:
     * intensity column = mass component (spectrum) with any column name
     * 'mz' and 'level' columns
    :param scan_id: Globally unique scan ID to assign to the mzXML <scan> structure
    :return: String containing the mzXML <scan> structure for the given mass component
    """
    spectrum_ms1 = spectrum[spectrum['level'] == 1].drop(columns='level')
    spectrum_ms2 = spectrum[spectrum['level'] == 2].drop(columns='level')
    component_idx = spectrum_ms1.drop(columns='mz').columns.values[0]

    if not spectrum_ms1.empty:
        max_intensity_idx = spectrum_ms1[component_idx].idxmax()
        pseudo_precursor = {'prec_mz': spectrum_ms1.loc[max_intensity_idx, 'mz'],
                            'prec_intensity': spectrum_ms1.loc[max_intensity_idx, component_idx]}
    else:
        pseudo_precursor = {'prec_mz': isolation_window_center + np.random.random() * 0.1,
                            'prec_intensity': 300 + np.random.random()}  # Arbitrary value based on some inspected scans

    spectrum_ms2 = spectrum_ms2.sort_values('mz')
    spectrum_ms2 = [spectrum_ms2['mz'].values,
                    spectrum_ms2.drop(columns='mz').values.T[0]]

    if len(spectrum_ms2) == 0:
        logger.warning('MS 2 data in component %d has no intensity values above threshold' % component_idx)
        return ''

    try:
        xml_content = spectrum_to_xml(spectrum=spectrum_ms2,
                                      level=2,
                                      scan_no=scan_id,
                                      pseudo_precursor=pseudo_precursor)

        return xml_content

    except ValueError as e:
        logger.warning(str(e))
        logger.warning(spectrum_ms2)
        return ''


def spectrum_to_xml(spectrum, level, scan_no, pseudo_precursor):
    """
    :param spectrum: List of 2 numpy arrays: mz values and intensities
    :param scan_no: Scan ID assigned to this spectrum (globally unique in end file).
    :param pseudo_precursor: Dictionary with "prec_intensity" and "prec_mz" of this spectrum
    :return: String with mzXML encoding of input spectrum
    """
    xml_content = '<scan num="{}" '.format(scan_no)
    xml_content += 'scanType="Full"\n'
    #xml_content += 'centroided="0"\n'
    xml_content += 'msLevel="{}"\n'.format(level)
    xml_content += 'peaksCount="{}"\n'.format(len(spectrum[0]))
    # xml_content += 'polarity="+"\n')
    # xml_content += 'retentionTime="PT{}S"\n'.format(rt))
    xml_content += 'lowMz="{}" '.format(np.min(spectrum[0]))
    xml_content += 'highMz="{}" '.format(np.max(spectrum[0]))
    xml_content += 'basePeakMz="{}"\n'.format(
        spectrum[0][np.argmax(spectrum[1])])
    xml_content += 'basePeakIntensity="{}"\n'.format(np.max(spectrum[1]))
    xml_content += 'totIonCurrent="{}"\n'.format(sum(spectrum[1]))
    xml_content += '>\n'
    if level == 2:
        prec_attr = '<precursorMz precursorIntensity="{0}">{1}</precursorMz>\n'
        xml_content += prec_attr.format(pseudo_precursor['prec_intensity'],
                                        pseudo_precursor['prec_mz'])
    xml_content += '<peaks compressionType="none"\n'
    xml_content += 'compressedLen="0"\n'
    xml_content += 'precision="32"\n'
    xml_content += 'byteOrder="network"\n'
    xml_content += 'contentType="m/z-int">'
    xml_content += as_base64_string(spectrum)
    xml_content += '</peaks>\n'
    xml_content += '</scan>\n'
    return xml_content


def record_scan_count(mzxml_filename):
    """
    Update the <msRun scanCount> attribute.
    Done at the end since the mzXML is created from multiple sources.
    """
    logger.info('Recording scan count...')
    nscans = 0
    with open(mzxml_filename, 'r') as xml_file:
        for line in xml_file:
            if '<scan' in line:
                nscans += 1

    xml_contents = []
    with open(mzxml_filename, 'r') as xml_file:
        for line in xml_file:
            if line.startswith('<msRun'):
                xml_contents.append(f'<msRun scanCount="{nscans}">\n')
            else:
                xml_contents.append(line)

    with open(mzxml_filename, 'w') as xml_file:
        xml_file.writelines(xml_contents)

    logger.info('...done.')


def index_mzxml_file(mzxml_filename):
    logger.info('Indexing file...')

    with pyteomics_mzxml.MzXML(mzxml_filename) as reader:
        offsets = reader._offset_index

    with open(mzxml_filename, 'r') as orig_file:
        lines = orig_file.readlines()[:-1]

    lines.append('<index name = "scan">\n')
    for idx, offset in offsets['scan'].index_sequence :
        lines.append(f"<offset id = \"{idx}\">{offset}</offset>\n")
    lines.append('</index>\n')

    with open(mzxml_filename, 'w') as indexed_file:
        indexed_file.writelines(lines)

    with open(mzxml_filename, 'rb') as indexed_file:
        for line in indexed_file:
            if line.decode('ascii').startswith('<index'):
                index_offset = indexed_file.tell() - len(line)
                break

    with open(mzxml_filename, 'a') as indexed_file:
        indexed_file.write(f'<indexOffset>{index_offset}</indexOffset>\n')
        indexed_file.write('</mzXML>\n')

    logger.info('...done.')


def convert_mzxml_to_mgf(input_mzxml, output_mgf):
    """
    http://www.matrixscience.com/help/data_file_help.html
    https://fiehnlab.ucdavis.edu/projects/lipidblast/mgf-files
    """
    spectra = []
    for spectrum in pyteomics_mzxml.read(input_mzxml):
        scan_no = spectrum['num']
        prec_mz = spectrum['precursorMz'][0]['precursorMz']
        prec_intensity = spectrum['precursorMz'][0]['precursorIntensity']
        spectrum['params'] = {'title': 'scan=' + scan_no,
                              'pepmass': (prec_mz, prec_intensity)}
        spectra.append(spectrum)
    mgf.write(spectra=spectra,
              output=output_mgf,
              fragment_format='%.{0}g %.{1}g'.format(MZ_DIGITS, I_DIGITS),
              use_numpy=True,
              file_mode='w')


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
        slice_df = slice_df.melt(var_name='mz', value_name='intensity')
        slice_df['cycle'] = slice_df.index
        # TODO: Spliting is slow
        mz_level = slice_df['mz'].str.split('_ms', expand=True)
        slice_df['level'] = mz_level[1].astype(np.uint8)
        slice_df['mz'] = mz_level[0]
        slice_df['sample_no'] = sample_no
        slice_df_per_sample.append(slice_df)
    slice_df_per_sample = pd.concat(slice_df_per_sample, ignore_index=True)
    return slice_df_per_sample


def as_dataframe_with_mz_index(mass_mode, mz_indices):
    """
    :param mass_mode: intensity vectors as columns (of equal length) named with component numbers
    :param mz_indices: list of m/z indices as 123.45_ms1
    """
    if not isinstance(mz_indices[0], str):
        raise Exception("m/z indices must be a string of form like 100.5_ms1")

    if not '_ms' in mz_indices[0]:
        raise Exception("m/z indices must be a string of form like 100.5_ms1")

    if not isinstance(mass_mode, pd.DataFrame):
        mass_mode = pd.DataFrame(mass_mode, dtype=np.float32)

    mz_indices = pd.Series(mz_indices, index=mass_mode.index)
    mz_series, level_series = mz_indices.str.split('_ms').str
    mass_mode = mass_mode.assign(mz=mz_series.astype(np.float32))
    mass_mode = mass_mode.assign(level=level_series.astype(np.int8))
    return mass_mode


def expand_list_of_mz_indices_to_df(mz_indices):
    """
    Take list of string mz_indices as e.g. ['1000.75_ms1', '100.5_ms2']
    and converts to pandas.DataFrame with columns ['mz', 'level']
    The index reflects the order in the input list.
    """
    mz_indices = pd.DataFrame(mz_indices)
    mz_indices = mz_indices.assign(mz=mz_indices[0].str.split('_', expand=True)[0].astype(np.float32))
    mz_indices = mz_indices.assign(level=mz_indices[0].str.split('ms', expand=True)[1].astype(np.uint8))
    mz_indices = mz_indices.drop(columns=0)
    return mz_indices


def mz_indices_for_level(mz_index_list, level):
    """Return list indices for the mz values at the given MS level"""
    lvl_marker = '_ms' + str(level)
    return [i for i, mz in enumerate(mz_index_list) if mz.endswith(lvl_marker)]


def as_base64_string(spectrum):
    peak_pairs = np.ravel(list(zip(spectrum[0], spectrum[1])))
    peak_pairs = np.array(peak_pairs, dtype=np.float32)
    # Data must be in "network" == big endian byte order
    # http://tools.proteomecenter.org/formats/mzXML/mzXML_xmlspy_docs.html
    # http://sashimi.sourceforge.net/schema_revision/mzXML_2.1/Doc/mzXML_2.1_tutorial.pdf
    if sys.byteorder == 'little':
        peak_pairs = peak_pairs.byteswap()
    spectrum_bytes = peak_pairs.tobytes()
    encoded_spectrum = base64.standard_b64encode(spectrum_bytes)
    return encoded_spectrum.decode('utf-8')


def get_swath_intervals_from_mzml(mzml_filename):
    logger.info('Creating swath intervals file from %s ...' % mzml_filename)
    msrun = _get_pymzml_reader(mzml_filename)

    if msrun:
        isol_windows = itree.IntervalTree()

        for spectrum in msrun:
            spectrum_index = int(spectrum._xmlTree.get('index'))
            level = spectrum.get('ms level', -1)
            if level is None:
                logger.warning('Spectrum %d has no MS level information. Skipping.' % spectrum_index)
                continue

            if (level == 2) and ("precursors" in spectrum):
                prec_mz = spectrum["precursors"][0]["mz"]
                prec_isolation_window_lower_offset = spectrum.get('MS:1000828', '')
                prec_isolation_window_upper_offset = spectrum.get('MS:1000829', '')
                if prec_isolation_window_lower_offset:
                    isol_windows.addi(prec_mz - prec_isolation_window_lower_offset,
                                      prec_mz + prec_isolation_window_upper_offset,
                                      spectrum_index)
                else:
                    msg = ('Spectrum index %d has '
                           'no precursor isolation window lower offset')
                    logger.warning(msg % spectrum_index)
        logger.info('...done.')
        return isol_windows

    else:
        errmsg = 'Could not read file ' + mzml_filename
        logger.error(errmsg)
        raise IOError(errmsg)


def _get_pymzml_reader(mzml_filename):
    return pymzml.run.Reader(mzml_filename,
                             extraAccessions=[('MS:1000016', ['value', 'unitName']),
                                              ('MS:1000828', ['value']),
                                              ('MS:1000829', ['value']),
                                              ('MS:1000040', ['value'])
                                              ])


def save_isolation_windows_to_csv(isol_windows, csv_filename):
    isol_windows = [tuple(itv)[0:2] for itv in isol_windows]
    isol_windows = sorted(set(isol_windows))
    with open(csv_filename, 'w') as csv_output:
        writer = csv.writer(csv_output)
        for window in isol_windows:
            writer.writerow(window)


def get_swath_intervals_from_adjusted_map_csv(filename):
    """Return the list of (adjusted) swath intervals in filename"""
    return get_swath_intervals_from_adjusted_map(pd.read_csv(filename))


def get_swath_intervals_from_adjusted_map(msmap):
    """Return the list of (adjusted) swath intervals, sorted by lower bound"""
    msmap = msmap.dropna(subset=['swath_lower_adjusted'])
    msmap = msmap.drop_duplicates('swath_lower_adjusted')

    swaths = zip(msmap['swath_lower_adjusted'], msmap['swath_upper_adjusted'])
    swaths = sorted(swaths, key=lambda x: x[0])
    return swaths


def get_filenames_with_extension(directory, extension):
    """
    Return list of relative paths and file names matching file extension in given
    directory.
    """
    extension = extension.split('.')[-1]  # Make sure to avoid duplicate dot
    filenames = []
    for filename in os.listdir(directory):
        if filename.endswith('.' + extension):
            filenames.append(os.path.join(directory, filename))
    return filenames


def get_swath_from_full_map_csv(swath_interval, filename):
    """
    Return only values within the swath from a full CSV map.

    :return: pandas.Dataframe with the corresponding entries in filename
    """
    msmap = pd.read_csv(filename)
    msmap = pd.concat([
        msmap[(msmap['level'] == 1) & (msmap['mz'].between(swath_interval[0],
                                                           swath_interval[1]))],
        msmap[(msmap['level'] == 2) & (msmap['prec_mz'].between(swath_interval[0],
                                                                swath_interval[1]))]
    ])
    return msmap


def get_swaths_from_csv_using_lookup_table(swath_no, swaths_dir, lookup_table):
    """
    Return only values within the swath from a full CSV map.

    :param lookup_table: pandas.Dataframe matching swath_no to swath_file
    :return: pandas.Dataframe with the corresponding entries in filename
    """
    swath_files = sorted(
        lookup_table[lookup_table['swath_no'] == swath_no]['swath_file'])
    swath_maps = [pd.read_csv(os.path.join(swaths_dir, f)) for f in swath_files]
    return swath_maps


def round_down_mz_and_sum_intensity(msmap_dataframe, decimals=2):
    m = msmap_dataframe.assign(mz=np.around(msmap_dataframe.mz,
                                            decimals=decimals))
    m = m.groupby(['level',
                   'rt',
                   'mz',
                   'prec_mz',
                   'prec_isolation_window_lower_offset',
                   'prec_isolation_window_upper_offset'],
                   as_index=False).sum()
    return m


def round_down_rt_and_avg_intensity(msmap_dataframe, decimals=4):
    m = msmap_dataframe.assign(rt=np.around(msmap_dataframe.rt,
                                            decimals=decimals))
    m = m.groupby(['level',
                   'rt',
                   'mz',
                   'prec_mz',
                   'prec_isolation_window_lower_offset',
                   'prec_isolation_window_upper_offset'],
                   as_index=False).mean()
    return m


def filter_out_MS2_entries_with_intensity_in_lower_percent(spectrum: pd.DataFrame,
                                                           intensity_cutoff_bin=0,
                                                           bins=100):
    """
    Discard MS 2 entries with intensity values that fall in the lower bins of the
    intensity distribution.

    Most values are expected to fall within this range
    and may be considered as background.

    :param spectrum:  mass_mode[[component_idx, 'mz', 'level']]
    :param bins: Number of bins to use when building intensity histogram
    :param intensity_cutoff_bin: The upper bound of this bin will be cutoff value
    :return: filtered pandas.DataFrame
    """
    spectrum_ms2 = spectrum[spectrum['level'] == 2]
    intensity_col = list(set(spectrum.columns.values) - set(['mz',  'level']))[0]

    intensity_hist = np.histogram(spectrum_ms2[intensity_col], bins=bins)
    intensity_cutoff = intensity_hist[1][intensity_cutoff_bin]

    filtered_spectrum = spectrum[
        (spectrum['level'] == 1) |
        ((spectrum['level'] == 2) & (spectrum[intensity_col] > intensity_cutoff))
    ]
    return filtered_spectrum


def undo_min_max_feature_scaling(data, min_val, max_val):
     return data * (max_val - min_val) + min_val


def get_filename_with_most_time_points(sample_file_list):
    with Pool(processes=cpu_count()) as pool:
        time_lengths = pool.map(count_time_points, sample_file_list)
    return sample_file_list[np.argmax(time_lengths)]


def count_time_points(sample_filename):
    time_points = []
    with open(sample_filename, 'r', newline='') as sample_csv:
        sample_reder = csv.reader(sample_csv, delimiter=',')
        for row in sample_reder:
            if sample_reder.line_num > 1:
                time_points.append(row[1])

    return len(set(time_points))


def all_fragments(peptide, charge, types=('b', 'y')):
    """
    Generate all possible m/z for fragments of given types and charge.
    From https://pyteomics.readthedocs.io/en/latest/examples/example_msms.html
    """
    for length in range(1, len(peptide)-1):
        for ion_type in types:
            if ion_type[0] in 'abc':
                yield mass.fast_mass2(peptide[:length], ion_type=ion_type, charge=charge)
            else:
                yield mass.fast_mass2(peptide[length:], ion_type=ion_type, charge=charge)


def mz_tol(mz, machine_ppm=40):
    return (mz * machine_ppm) / 1e6


def index_of_nearest_value(series, value):
    return (np.abs(series - value)).idxmin()


def isin_float(data, points_array, tolerance=1e-5):
    """
    Alternative to pandas.DataFrame.isin() with float tolerance.

    Return Boolean values for entries in data (pandas.Series or numpy.array)
    that are close to any element in points (numpy.array)

    Note:
        Relative tolerance = 0 for equivalence with (x +/- epsilon) formulation.
    """
    return np.isclose(data,
                      points_array[:, np.newaxis],
                      atol=tolerance,
                      rtol=0).any(axis=0)


def align_columns(dataframe_list):
    """
    Return list of dataframes that have the same column indices.
    Empty columns are added at a previously non-existing index.
    """
    aligned_maps = []
    common_col_names = common_columns(dataframe_list)

    # Inserting columns in dataframes is really slow (by design).
    # For each map, we preallocate a larger one with all columns,
    # then insert the map into the larger one.
    for i in range(len(dataframe_list)):
        full_map = pd.DataFrame(columns=common_col_names)  # "headers" only
        full_map = full_map.append(dataframe_list[i])

        # Sorting because of Pandas bug #4588 (appending reorders columns)
        sorted_columns = sorted(full_map.columns)
        full_map.reindex(sorted_columns, axis=1, copy=False)
        aligned_maps.append(full_map)

    return aligned_maps


def common_columns(dataframe_list):
    common_col_names = set()
    for i in range(len(dataframe_list)):
        common_col_names |= set(dataframe_list[i].columns.values)
    return sorted(common_col_names, key=mz_label_sort_value)


def mz_label_sort_value(col_label):
    """
    For sorting m/z labels according to numerical value and MS level.
    """
    val, level = col_label.split('_')
    if level == "ms2":
        return float(val) * 100
    else:
        return float(val)


def get_chromatograms_at_mz_and_level(map_list, mz, level, time_range=None):
    """
    :return: List of pandas.DataFrame with columns 'rt' and 'intensity'
    """
    chromatograms = []
    for i in range(len(map_list)):
        chrom = map_list[i][(map_list[i]['level'] == level) &
                            (np.isclose(map_list[i]['mz'], mz))].sort_values('rt')

        if time_range is not None:
            chrom = chrom[chrom['rt'].between(time_range[0], time_range[1])]

        chromatograms.append(chrom)
    return chromatograms
