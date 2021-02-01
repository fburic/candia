# CANDIA: Canonical Decomposition of Data-Independent-Acquired Spectra

![logo](candia_logo.png)

CANDIA is a GPU-powered unsupervised multiway factor analysis framework that deconvolves multispectral scans to individual analyte spectra, chromatographic profiles, and sample abundances, using the PARAFAC (or canonical decomposition) method. The deconvolved spectra can be annotated with traditional database search engines or used as a high-quality input for *de novo* sequencing methods.

> *Parallel Factor Analysis Enables Quantification and Identification of Highly Convolved Data-Independent-Acquired Protein Spectra*  
> Filip Buric, Jan Zrimec, Aleksej Zelezniak, *Cell Patterns* (2020); DOI: https://doi.org/10.1016/j.patter.2020.100137

## Note

This repository holds the CANDIA scripts that produced the results in the submitted article 
(in the revision with the tag ["submission"](https://github.com/fburic/candia/releases/tag/submission)).
Restructuring the pipeline and making it easy to use is a *work in progress*. 
The repository will be updated but the `submission` revision is available as a snapshot of the source code 
at the time of manuscript submission.


## Hardware Requirements

For the decomposition, **an NVIDIA GPU is required**. Models known to perform well: K80, GP100, and V100.
Minimum GPU RAM: 8 GB. Recommended: >= 16 GB.

The preprocessing and downstream stages should perform acceptably with as few as 8 CPUs at 2.5 GHz and 16 GB RAM.
Recommended: >= 16 CPUs at >= 3 GHz, and >= 32 GB RAM.


# User guide

- [Installation](#installation)
    + [Third-party software](#third-party-software)
- [Usage](#usage)


## Installation

The current CANDIA distribution is split in two components: 

* the scripts that make up the pipeline (clone this repo to fetch them)
* a Singularity container which includes all dependencies for running CANDIA (Python libraries and other necessary third-party software)

We recommend using the Singularity container provided on Singularity Hub (built with Singularity version 3.4.2). 
The container should be placed inside the cloned CANDIA repo.
To use the container:

1. [Install](https://sylabs.io/guides/3.0/user-guide/quick_start.html) Singularity 3.x on your system
2. Pull the container with the command `singularity pull candia.sif shub://fburic/candia:def` This will download the container `candia.sif` (4.4 GB) inside the current working directory. The main `candia` script assumes `candia.sif` is located in the CANDIA repo.

Alternatively, the container may be built from the supplied `candia.def` file 
(requires root permissions, see instructions [here](https://sylabs.io/guides/3.0/user-guide/build_a_container.html)).

The bulk of dependencies is managed through the conda packaging system. 
The container simply encapsulates such an enviroment to ensure portability, but 
a conda environment can also be [built](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file) 
from scratch on the user's system with the provided `candia_env.yaml` specification file.
Here, commands will be shown using the Singularity container.

CANDIA was developed and used on POSIX systems: Ubuntu 18.04 and CentOS 7.8 Linux, 
as well as macOS 10.13 (High Sierra), 10.14 (Mojave), 10.15 (Catalina).
By using the Singularity container, CANDIA should be runnable on any OS that supports Singularity.


### Third-party software

The Singularity container (and associated conda environment) currently includes:

* Crux 3.2 (package `crux-toolkit`, `bioconda` channel) 
* TPP 5.0.0 (package `tpp`, `bioconda` channel)
* msproteomicstools 0.11.0 (`bioconda` channel)

The software below are either not distributed through package management archives or 
those versions did not work with this setup.

#### MS-GF+

MS-GF+ comes in the form of a Java JAR package and it only needs to be unzipped into
a directory (e.g. `/home/$USER/software`).
Installation instructions: https://github.com/MSGFPlus/msgfplus

The environment variable `MSGF_JAR_PATH` needs to be set to inform the pipeline 
of the `.jar` location. Add it to your `.bashrc` or `.profile`. E.g.

```shell script
export MSGF_JAR_PATH="$HOME/software/MSGFPlus/MSGFPlus.jar"
```

#### DIA-NN

To perform quantification with the generated CANDIA library, we provide a wrapper script
for the Linux version of DIA-NN. The the user may install `diann-linux` from https://github.com/vdemichev/DiaNN 
(included in the "source code" package of each release).

CANDIA is known to work with DIA-NN `1.7.4`

The wrapper `scripts/quantification/diann.Snakefile` is simply a convenience script that supplies all necessary parameters,
and ensures resuming on error.
Once the library is created, the user may run DIA-NN independently of CANDIA or the Singularity workflow.


#### Mayu

Note: The `tpp` version `5.0.0-0` bioconda package that also includes Mayu may give missing error libraries.
The user can install a stand-alone version as described here.

Mayu is in the form of Perl scripts, so it only needs to be unzipped to be used.
http://proteomics.ethz.ch/muellelu/web/LukasReiter/Mayu/

It was designed to be executed from its installation directory 
(e.g. `/home/$USER/software/Mayu`) so the environment variable `MAYU_STANDALONE_PATH` 
needs to be set up to point at this location, so the pipeline may run it from anywhere.

Add it to your `.bashrc` or `.profile`. E.g.

```shell script
export MAYU_STANDALONE_PATH="$HOME/software/Mayu"
```



## Usage

The pipeline currently consists of a collection of scripts for the different stages of processing.
The shell script `candia` is provided to run all pipeline steps from preprocessing 
up to and including quantification library creation.
However, *this script is still a work in progress* and is currently mainly used to test that the pipeline
can execute on the user's system. 
This can be done by running the following from CANDIA's top level directory:

```shell script
./candia test/test_experiment/config/candia.yaml
``` 

The general syntax of this script is `candia EXPERIMENT_CONFIG_YAML`

**Note** The test experiment will only run properly up to and including the decomposition step.
As this is the most complex part of the pipeline, this should guarantee that your installation
works for real data. A more complete test suite (more realistic data and unit tests) is planned.

Currently, it is **recommended to run each stage at a time**, using the corresponding scripts. 
Step-by-step instruction are listed here using the toy data files provided in this repo. 
These data files are primarily meant to test the pipeline is working properly.

Expected running times are given as a very rough estimate on 9 real DIA scan files,
to indicate the time scale of the steps.
These times will of course vary depending on the available resources.
Actual running times for the test experiment should be exceedingly short.

We recommend running the pipeline on a workstation with multiple cores or 
a high-performance computing (HPC) environment. Currently, the HPC manager `slurm`
is supported. A conda environment is also assumed to be present on the cluster.
A Singularity recipe to encapsulate the environment is provided. 

An NVIDIA GPU card is required for the decomposition.
Partial support for CPU-only execution is implemented but the performance becomes infeasible.


### 0. Configure CANDIA execution

The execution of the pipeline is **configured through a YAML file**
(e.g. `test_experiment/config/candia.yaml`). 
This configuration file specifies the location of input and intermediate files,
as well as algorithm parameters. The paths are interpreted as relative to the 
experiment directory `root_dir`.

#### Good practice

It is recommended to create **separate configuration files** for each parameter variation to be included in final results,
rather than changing a single configuration file.
Different configuration files may also be created for each stage of the pipeline.
The important thing is to supply the relevant parameters to each script.

Please read through the next steps to see which parameters are relevant at each step.
A full configuration file is provided in the test experiment.
More information may be found in the CANDIA article Supplemental Information.

The scripts are **expected to be run from inside the CANDIA repo** 
(Thus first `cd candia` after clonning).


### 1. Convert DIA scan files from mzML to CSV

> Expected running time: 3-5 minutes per file

The basename of the input mzML files will be kept downstream, with only the extension changed.

Relevant pipeline config values:

* `samples_mzml` - location of input mzML files
* `samples_csv` - location of output CSV files
* `swath_windows` - file listing precursor isolation windows 
* `min_scan_intensity` - all values below this are dropped

Note: To set Snakemake to use `N` cores, pass it the `--jobs N` argument.
It should automatically use all available cores. 

Commands:

```shell script
configfile='test/test_experiment/config/candia.yaml'

singularity exec candia.sif \
    snakemake -s scripts/util/mzml2csv.Snakefile --configfile ${configfile}
```


### 2. Adjust precursor isolation windows

We cannot have overlapping isolation windows, so the bounds are adjusted.

> Expected running time: 2 min per file

Relevant pipeline config values:

* `samples_csv` - location of CSV scan files
* `samples_adjusted_swaths` - location of output adjusted CSV files
* `swath_windows_adjusted` - file listing adjusted precursor isolation windows

Command:

```shell script
singularity exec candia.sif \
    snakemake --forceall -p -s scripts/util/adjust_swaths.Snakefile --configfile  ${configfile}
```

Save the adjusted intervals for downstream tasks:

```shell script
EXP_DIR=$(grep "root_dir:" ${configfile} | cut -d " " -f 2 | tr -d \")
cp $(find ${EXP_DIR} -name "*.intervals" | head -n 1) \
   ${EXP_DIR}/$(grep "swath_windows_adjusted:" ${configfile} | cut -d " " -f 2 | tr -d \")
```


### 3. Split samples into slices

The input scans are partitioned into swaths and RT windows of width `window_size_sec`.
Note that this will create many small CSV files in `slices_location/swath=value/rt_window=value` subdirectories. 

> Expected running time: 30 min

Relevant pipeline config values:

* `samples_adjusted_swaths` - location of output adjusted CSV files
* `window_size_sec` - the width of the RT windows
* `slices_location` - location of output slices

```shell script
singularity exec candia.sif \
    python scripts/util/split_csv_maps_to_slices.py --config ${configfile}
```


### 4. Generate tensor files for all slices

The slices are converted into NumPy (swath, rt, sample) tensors stored as `.npy` files. 
The Snakefile will start the conversions for each slice in parallel and each such job may
be submitted independently as a HPC cluster job.

> Expected running time: 10 - 40 min (depending on available resources)

Relevant pipeline config values:

* `mass_tol_ppm` - The mass tolerance (in PPM) of the scan acquisition

The m/z precision is 10 decimals and 
m/z partitions with less than 5 time points in any sample are filtered out.

```shell script
singularity exec candia.sif \
    snakemake --jobs 4 -s scripts/util/generate_slice_tensors.Snakefile \
    -R generate_slice_tensor --forceall --rerun-incomplete --configfile ${configfile} 
```

On a slurm-managed cluster, this would be:

```shell script
NJOBS=200
ACCOUNT_NUMBER=ABC123
singularity exec candia.sif \
    snakemake --jobs ${NJOBS} -s scripts/util/generate_slice_tensors.Snakefile \
    -R generate_slice_tensor --forceall --rerun-incomplete --configfile ${configfile} \
    --cluster "sbatch -A ${ACCOUNT_NUMBER} -t 06:00:00 --ntasks 6"
```


### 5. Run PARAFAC decomposition

A decomposition is performed on each slice tensor, for all number of components `F`
in the configured range. Multiple tensors are decomposed in parallel on each available GPU card.
If multiple GPUs ara available, the set of input tensors is partitioned evenly between the cards.

Note that while the scripts below require a GPU, the decomposition script itself 
`scripts/parafac/decompose_parafac.py` may be run on CPUs if executed without the `--use_gpu` flag. 

Relevant pipeline config values:

* `parafac_*` - decomposition parameters
* `parafac_min_comp` - Lower bound of the number of components `F` to decompose for
* `parafac_max_comp` - Upper bound of the number of components `F` to decompose for

> Expected running time: 6 - 12 hours (depending on the GPU model)

#### Workstation

The workstation command has the syntax:

`scripts/parafac/decompose_workstation.sh  EXPERIMENT_CONFIG_FILE N_PARALLEL_DECOMP_PER_GPU [SNAKEMAKE_OPTIONS]`

Example running with Singularity (with 2 parallel decompositions per GPU):

```shell script
scripts/parafac/decompose_workstation.sh ${configfile} 2
```

#### HPC cluster

The cluster command has the syntax:

`scripts/parafac/decompose_cluster.sh ${snake_configfile} N_PARALLEL_DECOMP_PER_GPU`

And is designed to be run as an array job on slurm.
Example (with 6 parallel decompositions per GPU):

```shell script
sbatch --array=0-1 --gres=gpu:1 --time 24:00:00 \
    scripts/parafac/decompose_cluster.sh "${configfile}" 6
```

If Singularity is not available, you can use the
`decomopse_cluster_no_singularity.sh` script. 


### 6. Index all PARAFAC models and components

All PARAFAC models and components are indexed with a unique ID, as support for
downstream tasks. These IDs, along with model filenames, are saved in two database
(or "index") files, in [Apache Feather](https://arrow.apache.org/docs/python/feather.html) format.
 
> Expected running time: a few sec
 
Relevant pipeline config values:

* `model_index` - Model ID index filename (Apache Feather format)
* `spectrum_index` - PARAFAC component ID index filename (Apache Feather format)
 
```shell script
singularity exec candia.sif \
    python scripts/parafac/models.py -c ${configfile}
```


### 7. Select Best Models

Measure the unimodality of all PARAFAC components (for all models),
then create a list of the best PARAFAC models, according to the unimodality criterion. 

> Expected running time: 10 sec

Relevant pipeline config values:

* `avg_peak_fwhm_sec` - The expected peak FWHM, to inform the peak finding procedure.
* `window_size_sec` - Width of the RT windows (in seconds).
* `time_modes_values` - Tabular file in Feather format that collects peak counts for all time modes (in all components), 
  for all models generated by CANDIA across all slices.
* `best_models` - CSV file where the best model parameter values are stored

First, collect the time modes and measure unimodality fractions.

```shell script
singularity exec candia.sif \
    python scripts/parafac/collect_time_mode_values.py -c ${configfile}
```

Then, select models:

```shell script
singularity exec candia.sif \
    Rscript scripts/parafac/select_best_models.R -c ${configfile}
```


### 8. Identify proteins using Crux or MS-GF+

Configure which tool to use by setting the `analysis_pipeline` 
variable to either `"crux"` or `"msgf+"` in the configuration file.
For example, in the provided test experiment, this is set to `"crux"`

A separate decoy sequence database must be provided besides the targets.

> Expected running time: 30 min with Crux, 1 h with MS-GF+ (depends on the set number of modifications)

Relevant pipeline config values:

* `best_models_mzxml` - MzXML file where the spectra from the best models are concatenated
* `best_models_crux_out_dir` - The directory for Crux results (Crux outputs multiple files)
* `database` - FASTA protein sequence database 
* `decoy_database` - FASTA decoy sequence database 
* `decoy_prefix` - Decoy sequence ID prefix, which should include separator (e.g. `decoy_`)
* `msgf_modifications` - MS-GF+ modifications configuration file, if any
* `msgf_threads` - The number of MS-GF+ computation threads to use

Crux version 3.2 is included in the Singularity container (and associated conda environment)

```shell script
singularity exec candia.sif \
    python scripts/identification/id_models_concat.py -c ${configfile}
```

If running CANDIA with MS-GF+, the container may be executed with a supplied path environment variable
specifying the location of the MS-GF+ program.

```shell script
MSGF_JAR_PATH="$HOME/software/MSGFPlus/MSGFPlus.jar" singularity exec candia.sif \
    python scripts/identification/id_models_concat.py -c ${configfile}
```


### 9. Build library

> Expected running time: 2-5 min

The [Schubert et. al (2015)](https://www.nature.com/articles/nprot.2015.015) 
protocol has been implemented and can be run as below.
For a more detailed description, please see the Supplemental Information of the CANDIA paper. 
The relevant pipeline config values are also described there.

A mixed target-decoy is required. This is specified through the `mixed_database` parameter

```shell script
singularity exec candia.sif \
    snakemake -s scripts/quantification/build_library.Snakefile --configfile ${configfile}
```

**Note** This protocol will fail for the supplied toy test data.


### 10. Quantify proteins with DIA-NN

To quantify proteins using the generated CANDIA library, a wrapper script is provided for running DIA-NN.
The path to the `diann-linux` binary should be supplied to the Singularity container (see command below).

Relevant pipeline config values:

* `samples_mzml` - The DIA scan files
* `quant_library` - The library file created by CANDIA
* `database` - The target protein sequence FASTA database
* `diann_out_dir` - Directory for DIA-NN output
* `diann_report` - The filename of the TSV report output by CANDIA

> Expected running time: 2-5 min per scan file

```shell script
SINGULARITYENV_PREPEND_PATH=$HOME/software/diann singularity exec candia.sif \
    snakemake -p -s scripts/quantification/diann.Snakefile --configfile ${configfile}
```


### 11. De novo sequencing with Novor and DeepNovo

Configure which tool to use through the configuration file

> Expected running time: 25 min with Novor, 10 min with DeepNovo (depends on GPU parameters) 

TODO: expand, clarify

```bash
snakemake -s scripts/denovo/sequence_best_models.Snakefile --configfile ${configfile}
```
