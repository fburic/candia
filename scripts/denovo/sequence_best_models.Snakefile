import inspect
import os
import sys
from pathlib import Path

_this_filename = inspect.getframeinfo(inspect.currentframe()).filename
_this_path = Path(_this_filename).parent.resolve()
sys.path.append(str(_this_path.parent))

import seqproc


blastp_columns = '\t'.join(['qseqid','sseqid','pident','length',
                            'mismatch','gapopen','qstart','qend',
                            'sstart','send','evalue','bitscore'])


rule all:
    input:
        config['sequencer_output_file']


rule mzxml2mgf:
    input:
        config['best_models_mzxml']
    output:
        config['best_models_mzxml'] + '.mgf'
    log:
        str(Path(config['log_dir']) / 'mzxml2mgf.log')
    shell:
        """
        FileConverter -in {input} -out {output} > {log} 2>&1

        python scripts/util/adjust_mgf_files.py -i {output} -o {output}.tmp \
        --sequencer {config[sequencer]}> {log} 2>&1 \
        && mv {output}.tmp {output}
        """


rule sequencer:
    input:
        rules.mzxml2mgf.output
    output:
        config['sequencer_output_file']
    log:
        str(Path(config['log_dir']) / 'denovo_sequencer.log')
    run:
        if config['sequencer'] == 'deepnovo':
            input_path = os.path.abspath(str(input))
            output_dir = os.path.dirname(os.path.abspath(str(input)))
            output_file = os.path.basename(str(output))
            shell('set +u ; source deactivate && source activate deepnovo ; export DEEPNOVO_INPUT={input_path} '
                  '&& pushd ${{USER}}/software/DeepNovo '
                  '&& python deepnovo_main.py --train_dir train.example --decode --beam_search --beam_size 5 '
                  '&& cp train.example/decode_output.tab {output_dir}/{output_file} '
                  '&& popd '
                  '&& source deactivate' )

        elif config['sequencer'] == 'novor':
            shell('novor.sh -f -p {config[novor_param_file]} -o {output} {input} > {log} 2>&1')

        else:
            raise NotImplementedError
