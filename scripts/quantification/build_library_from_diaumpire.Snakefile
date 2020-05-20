import os
from pathlib import Path
import yaml


# Set up us some paths
diaumpire_out_dir = Path(config['diaumpire_out_dir'])
lib_dir = diaumpire_out_dir / 'library'
lib_basename = lib_dir / 'diaumpire_library.xml'
mayu_path = os.getenv('MAYU_STANDALONE_PATH', default='')


# Disable PTMs for MS-GF+
config['msgf_modifications'] = ''


# Pass the YAML config to the rules that need it by writing it to a temp file
tmp_configfile_copy = f"{os.getenv('TMPDIR')}/config.yaml"
with open(tmp_configfile_copy, 'w') as configfile_copy:
    yaml.dump(config, configfile_copy, default_flow_style = False)


input_files = {
    fname.stem : str(fname) for fname in diaumpire_out_dir.glob('*.mgf')
}
def umpire_sources(wildcards):
    return input_files[wildcards.feature_file]


rule all:
    input:
        config["diaumpire_library"]


rule mgf2mzxml:
    input:
        str(diaumpire_out_dir / '{feature_file}.mgf')
    output:
        str(lib_dir / '{feature_file}.mzXML')
    shell:
        """
        FileConverter -in {input} -out {output}
        """


rule msgf:
    """
    Use MS-GF+ to get PSMs from DIA-Umpire pseudo-spectra (feature files)
    """
    input:
        rules.mgf2mzxml.output
    output:
        str(lib_dir / '{feature_file}.mzid')
    shell:
        """
        python scripts/util/wrappers.py --process msgf_mzid \
        -i {input} -o {output} \
        --config {tmp_configfile_copy}
        """


rule mzid2pepxml:
    input:
        rules.msgf.output
    output:
        touch(str(lib_dir / '{feature_file}.pep.xml'))
    shell:
        """
        idconvert {input} --pepXML --ext .pep.xml -o {lib_dir}
        """


rule peptide_prophet:
    """
    Use PeptideProphet to assign significance to mixed target-decoy Comet PSMs
    Run ProteinProphet afterwards
    """
    input:
        expand(rules.mzid2pepxml.output, feature_file = input_files.keys())
    output:
        lib_dir / 'interact.msgf.pep.xml'
    shell:
        """
        xinteract -OARPdp -d{config[decoy_prefix]} -N{output} {input}
        """


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
        mzid2pepxml_out = expand(rules.mzid2pepxml.output, feature_file = input_files.keys()),
        pepprophet_out = rules.peptide_prophet.output
    output:
        str(lib_dir / 'mayu.DONE')
    shell:
        """
        perl -I{mayu_path}/lib {mayu_path}/Mayu.pl \
        -verbose -A {input.pepprophet_out} -C {config[mixed_database]}  -E {config[decoy_prefix]} \
        -G {config[quant_lib_mayu_fdr]} \
        -P mFDR={config[quant_lib_mayu_fdr]}:t \
        -H 51 \
        -M {lib_dir}/$(basename {input.mzid2pepxml_out} .pep.xml) \
        && touch {output}
        """


rule spectrast:
    """
    Create library from significant PSMs and scan file.
    Using flag file for rule since one can only give spectrast a basename for its output files.
    http://tools.proteomecenter.org/wiki/index.php?title=Software:SpectraST#SpectraST_Options

    CUTOFF = The minimum iProphet probability at which the protein FDR is below threshold
    [Schubert et. al.]
    """
    input:
        rules.mayu.output,
        pepprophet_out = rules.peptide_prophet.output
    output:
        touch(str(lib_basename) + '.sptxt')
    shell:
        """
        set +o pipefail

        MAYU_OUT="{lib_dir}/_psm_mFDR{config[quant_lib_mayu_fdr]}_t_1.08.csv"
        CUTOFF=$(cat ${{MAYU_OUT}} | tail -n+2 | tr ',' '\\t' | cut -f5 | sort -g | head -n 1)

        spectrast -cN{lib_basename} \
        -cICID-QTOF \
        -cP${{CUTOFF}} \
        -c_RDY \
        -co \
        {input.pepprophet_out}
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
        spectrast_lib = rules.spectrast.output
    output:
        config["diaumpire_library"]
    shell:
        """
        set +o pipefail

        spectrast2tsv.py \
        -l {config[lower_mz_frag]},{config[upper_mz_frag]} \
        -s b,y \
        -x 2,3 \
        -o 4 -n 6 \
        -p {config[quant_library_spectrast_max_frag_annot_err]} \
        -d \
        -w <(cat {config[swath_windows]} | sed 's/,/\\t/g') \
        -k openswath -a {output} {input.spectrast_lib} || true
        """
