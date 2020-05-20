import random
import logging
import re

import pandas as pd
from pyopenms import *


logger = logging.getLogger(__name__)


def convert_sequencing_results_to_fasta(sequencing_results_filename: str,
                                        output_fasta_filename: str,
                                        sequencer: str,
                                        id_column: str = 'scan'):
    seq_results = read_sequencing_results(sequencing_results_filename, sequencer)
    save_sequence_dataframe_as_fasta(seq_results,
                                     sequence_column='sequence',
                                     header_column=id_column,
                                     fasta_filename=output_fasta_filename)


def read_sequencing_results(sequencing_results_filenames: str, sequencer: str) -> pd.DataFrame:
    """
    Read output file from the configured sequencer and
    return standard pandas.DataFrame with columns:   scan, sequence, score
    """
    if sequencer == 'deepnovo':
        seq_res = pd.read_csv(os.path.expandvars(sequencing_results_filenames), sep='\t')
        seq_res = seq_res[['scan', 'output_seq', 'output_score']]
        seq_res = seq_res[seq_res.output_seq != 'nan']
        # Fix seq strings
        seq_res['output_seq'] = seq_res['output_seq'].str.replace(',', '')
        # Use standard names
        seq_res = seq_res.rename(columns={'output_seq': 'sequence',
                                          'output_score': 'score'})
        # Remove modifications
        seq_res['sequence'] = seq_res['sequence'].str.replace('mod', '')

        return seq_res

    elif sequencer == 'novor':
        novor_columns = ['id', 'scanNum', 'RT', 'mz(data)', 'z',
                         'pepMass(denovo)', 'err(data-denovo)',
                         'ppm(1e6*err/(mz*z))',
                         'score', 'peptide', 'aaScore']
        novor_res = pd.read_csv(sequencing_results_filenames, skipinitialspace=True, comment='#',
                                names=novor_columns)
        # Use standard names
        novor_res = novor_res.rename(columns={'peptide': 'sequence',
                                              'scanNum': 'scan'})
        # Remove PTMs
        novor_res['sequence'] = novor_res['sequence'].str.replace(r'\(.+\)', '')

        return novor_res

    else:
        raise NotImplementedError


def get_evalue_threshold(sequencing_results_filename: str,
                         sequencer: str,
                         database_filename: str,
                         pvalue_threshold=0.05) -> int:
    """
    Get an E-value cutoff depending on the database size, median sequence length,
    and which corresponds to a given p-value (= 0.05, by default).
    """
    sequencing_results = read_sequencing_results(sequencing_results_filename, sequencer)
    median_seq_len = sequencing_results['sequence'].str.len().median()

    with open(database_filename, 'r') as db_file:
        database = db_file.read()
    num_db_sequences = sum([1 for _ in re.finditer('\n>', database)]) + 1

    evalue_threshold = pvalue_threshold * median_seq_len * num_db_sequences
    return int(np.ceil(evalue_threshold))


def save_sequence_dataframe_as_fasta(sequences: pd.DataFrame,
                                     sequence_column: str,
                                     header_column: str,
                                     fasta_filename: str):
    """
    Note: The file is first built in memory, then written, since I/O is a bottleneck.
    :param sequence_column: Which column holds the sequences
    :param header_column: Which column to use as FASTA headers
    """
    fasta_entries = []
    for header, sequence in sequences[[header_column, sequence_column]].values:
        fasta_entries.append(f'>{header}\n{sequence}\n')
    with open(fasta_filename, 'w') as output_fasta:
        output_fasta.writelines(fasta_entries)


def parse_blastp_hits(blast_reults_filename: str) -> pd.DataFrame:
    """
    Assumes blastp output is in format 6
    """
    columns = ['qseqid', 'sseqid', 'pident', 'length',
               'mismatch', 'gapopen', 'qstart', 'qend',
               'sstart', 'send', 'evalue', 'bitscore']
    blastp_res = pd.read_csv(blast_reults_filename,
                             sep='\t', header=None, names=columns)
    return blastp_res


def mutate_seq(sequence: str, mutation='random', percent=0.5) -> str:
    """
    Shuffle amino acids in each fragment produced by digestion.
    mutation can be 'rotate', 'random_end', 'random_center', random (default)
    """
    sequence = sequence.upper()
    sequence = AASequence.fromString(sequence)
    dig = ProteaseDigestion()
    frags = []
    dig.digest(sequence, frags)

    try:
        mutated_fragments = []
        for frag in frags:
            frag = frag.toString().decode('utf-8')
            if mutation == 'rotate':
                mutated_fragments.append(frag[:-5] + frag[-3:] + frag[-5:-3])

            elif mutation == 'random_center':
                mid_start = len(frag) // 2 - 3
                mid_end = len(frag) // 2 + 3
                mfrag = (frag[:mid_start]
                         + ''.join(random.sample(frag[mid_start:mid_end], k=len(frag[mid_start:mid_end])))
                         + frag[mid_end:])
                mutated_fragments.append(mfrag)

            elif mutation == 'random':
                positions = random.sample(range(len(frag) - 1), k=int(percent * (len(frag) - 1)))
                new_positions = np.random.permutation(positions).tolist()
                while len(positions) >= 2 and positions == new_positions:
                    new_positions = np.random.permutation(positions).tolist()
                mfrag = list(frag)
                for pos, new_pos in zip(positions, new_positions):
                    mfrag[new_pos] = frag[pos]
                mutated_fragments.append(''.join(mfrag))

            else:
                mutated_fragments.append(frag[:-5]
                                         + ''.join(random.sample(frag[-5:], k=len(frag[-5:]))))
        return ''.join(mutated_fragments)

    except IndexError as e:
        return sequence
