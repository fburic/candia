import os
from pathlib import Path

lib_dir = Path(config["root_dir"]) / config['quant_library_dir']
lib_basename = lib_dir / 'best_models_library.xml'
mayu_path = os.getenv('MAYU_STANDALONE_PATH', default='')


rule all:
    input:
        Path(config["root_dir"]) / config["quant_library"]


rule comet:
    """Run Comet on a mixed target-decoy library"""
    input:
        Path(config["root_dir"]) / config['best_models_mzxml']
    output:
        lib_dir / 'comet.target.pep.xml'
    shell:
        """
        crux comet --peptide_mass_units 2 --peptide_mass_tolerance {config[mass_tol_ppm]} \
         --overwrite T --output-dir {lib_dir} {input} {config[root_dir]}/{config[mixed_database]}
        """


rule peptide_prophet:
    """Use PeptideProphet to assign significance to mixed target-decoy Comet PSMs"""
    input:
        rules.comet.output
    output:
        touch(lib_dir / 'peptide_prophet.DONE')
    shell:
        "PeptideProphetParser {input} DECOY={config[decoy_prefix]}"


rule mayu:
    """
    Note: Need to use a standalone Mayu version. The TPP version is missing libraries.
    Flags:
    -G = maximal PSM FDR
    -H = number of analysis steps
    -P = output filtered ids
    -M = use this as file name base
    """
    input:
        rules.comet.output, rules.peptide_prophet.output
    output:
        touch(lib_dir / 'mayu.DONE')
    shell:
        """
        perl -I{mayu_path}/lib {mayu_path}/Mayu.pl \
        -verbose -A {input} -C {config[root_dir]}/{config[mixed_database]}  -E {config[decoy_prefix]} \
        -G {config[quant_lib_mayu_fdr]} \
        -P mFDR={config[quant_lib_mayu_fdr]}:t \
        -H 51 \
        -M {lib_dir}/$(basename {rules.comet.output} .pep.xml)
        """


rule make_scan_available:
    """
    SpectraST expects the scan file to have the 
    same basename as the PSM pepXML and be in the same dir
    """
    input:
        Path(config["root_dir"]) / config['best_models_mzxml']
    output:
        lib_dir / 'comet.mzXML'
    shell:
        "ln -s $(pwd)/{input} {output}"


rule spectrast:
    """
    Create library from significant PSMs and scan file.
    Using flag file for rule since one can only give spectrast a basename for its output files.
    http://tools.proteomecenter.org/wiki/index.php?title=Software:SpectraST#SpectraST_Options

    CUTOFF = The minimum iProphet probability at which the protein FDR is below threshold
    [Schubert et. al.]
    """
    input:
        rules.peptide_prophet.output,
        rules.make_scan_available.output,
        rules.mayu.output,
        rules.comet.output
    output:
        touch(str(lib_basename) + '.sptxt')
    shell:
        """
        set +o pipefail
        MAYU_OUT="{lib_dir}/$(basename {rules.comet.output} .pep.xml)_psm_mFDR{config[quant_lib_mayu_fdr]}_t_1.08.csv"
        CUTOFF=$(cat ${{MAYU_OUT}} | tail -n+2 | tr ',' '\\t' | cut -f5 | sort -g | head -n 1)

        spectrast -cN{lib_basename} \
        -cICID-QTOF \
        -cP${{CUTOFF}} \
        -c_RDY \
        -co \
        {rules.comet.output}
        """


rule spectrast2openswath:
    """
    -l = Lower and upper m/z limits of fragment ions
    -s = Ion types to consider
    -x = Charges to consider
    -n = Max number of reported ions per peptide/z. Default: 20
    -o = Min number of reported ions per peptide/z. Default: 3
    -p = Maximum error allowed at the annotation of a fragment ion. Default: 0.05
    -d = Remove duplicate masses from labeling
    -e = Use theoretical mass
    """
    input:
        spectrast_lib=rules.spectrast.output
    output:
        Path(config["root_dir"]) / config["quant_library"]
    shell:
        """
        spectrast2tsv.py \
        -l {config[lower_mz_frag]},{config[upper_mz_frag]} \
        -s b,y \
        -x 2,3 \
        -o 4 -n 6 \
        -p {config[quant_library_spectrast_max_frag_annot_err]} \
        -d \
        -w <(cat {config[root_dir]}/{config[swath_windows]} | sed 's/,/\t/g') \
        -k openswath -a {output} {input.spectrast_lib}
        """



