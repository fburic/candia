import inspect
import os
from multiprocessing import Pool, cpu_count
import sys
from pathlib import Path

import pandas as pd

_this_filename = inspect.getframeinfo(inspect.currentframe()).filename
_this_path = Path(_this_filename).parent.resolve()
sys.path.append(str(_this_path.parent))

from denovo import seqproc


blastp_columns = '\t'.join(['qseqid','sseqid','pident','length',
                            'mismatch','gapopen','qstart','qend',
                            'sstart','send','evalue','bitscore'])

# Set up us some paths
diaumpire_out_dir = Path(config['diaumpire_out_dir'])
out_dir = diaumpire_out_dir / (config['sequencer'] + '_out')

input_files = {
    fname.stem : str(fname) for fname in diaumpire_out_dir.glob('*_Q1.mgf')
}
def umpire_sources(wildcards):
    return input_files[wildcards.feature_file]


def load_results(filename):
    seq_res = seqproc.read_sequencing_results(filename, config['sequencer'])
    return seq_res.assign(filename=filename)

rule all:
    input:
        out_dir / 'denovo_diaumpire.csv'


rule sequencer:
    input:
        diaumpire_out_dir / '{feature_file}.mgf'
    output:
        touch(Path(config['tmp_dir']) / '{feature_file}_denovo.tsv')
    log:
        str(Path(config['log_dir']) / '{feature_file}_denovo.log')
    run:
        if config['sequencer'] == 'deepnovo':
            input_path = os.path.abspath(str(input))
            shell('set +u ; source deactivate && source activate deepnovo ; export DEEPNOVO_INPUT={input_path} '
                  '&& pushd /home/fburic/software/DeepNovo '
                  '&& echo > train.example/decode_output.tab '
                  '&& python deepnovo_main.py --train_dir train.example --decode --beam_search --beam_size 5 ; '
                  'cp train.example/decode_output.tab {output} ; '
                  'popd ; '
                  'source deactivate')

        elif config['sequencer'] == 'novor':
            shell('novor.sh -f -p {config[novor_param_file]} -o {output} {input} > {log} 2>&1')

        else:
            raise NotImplementedError


rule collate_sequences:
    input:
        expand(rules.sequencer.output, feature_file = input_files.keys())
    output:
        out_dir / 'denovo_diaumpire.csv'
    run:
        seq_results = []
        for infile in input:
            try:
                res = load_results(infile)
                seq_results.append(res)
            except Exception as e:
                print('File could not be read: ', str(infile))
        seq_results = pd.concat(seq_results, ignore_index=True)
        seq_results.to_csv(str(output), index=False)

