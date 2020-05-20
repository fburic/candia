from pathlib import Path

diann_input_files = Path(config['samples_mzml']).glob('*.mzML')
diann_input_string = ' '.join([f'--f "{fname}"' for fname in diann_input_files])
diann_out_dir = Path(config['diann_out_dir'])

diann_output_target = {
    'library-free': diann_out_dir / config['diann_library'],
    'normal': diann_out_dir / config['diann_report']
    }

run_type = config.get('runtype', 'normal')
    

rule all:
    input:
        diann_output_target[run_type]


rule diann_library_free:
    input:
        diann_input_files
    output:
        report = str((diann_out_dir / config['diann_report']).with_suffix('')) + '_libfree.tsv',
        library = diann_out_dir / config['diann_library'],
        gene_stats = (diann_out_dir / config['diann_report']).with_suffix('.genes.tsv')
    shell:
        """
        diann-linux {diann_input_string} \
        --out "{output.report}" \
        --out-gene "{output.gene_stats}" \
        --out-lib "{output.library}" \
        --fasta "{config[database]}" \
        --learn-lib "{config[diann_train_lib]}" \
        --gen-spec-lib \
        --lib "" \
        --qvalue 1 --fasta-search --min-fr-mz 200 --max-fr-mz 1800 --met-excision --cut-after KR --missed-cleavages 1 --min-pep-len 7 --max-pep-len 30 --min-pr-mz 300 --max-pr-mz 1800 --unimod4 --no-quant-files \
        --threads {config[diann_threads]} --verbose 3
        """


rule diann_quant:
    input:
        scans = diann_input_files,
        library = config['quant_library']
    output:
        report = diann_out_dir / config['diann_report'],
        gene_stats = (diann_out_dir / config['diann_report']).with_suffix('.genes.tsv')
    shell:
        """
        diann-linux {diann_input_string} \
        --lib "{input.library}" \
        --fasta "{config[database]}" \
        --out "{output.report}" \
        --out-gene "{output.gene_stats}" \
        --qvalue 1 --met-excision --no-quant-files \
        --threads {config[diann_threads]} --verbose 3
        """



