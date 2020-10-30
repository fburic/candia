# CANDIA: Canonical Decomposition of Data-Independent-Acquired Spectra

![logo](candia_logo.png)

CANDIA is a GPU-powered unsupervised multiway factor analysis framework that deconvolves multispectral scans to individual analyte spectra, chromatographic profiles, and sample abundances, using the PARAFAC (or canonical decomposition) method. The deconvolved spectra can be annotated with traditional database search engines or used as a high-quality input for *de novo* sequencing methods.

*Parallel factor analysis enables quantification and identification of highly-convolved data independent-acquired protein spectra*, Filip Buric, Jan Zrimec, Aleksej Zelezniak, bioRxiv 2020.04.21.052654; doi: https://doi.org/10.1101/2020.04.21.052654

## Note

This repository holds the CANDIA scripts that produced the results in the bioRxiv preprint. Restructuring the pipeline and making it easy to use is a work in progress. The repository will be updated but this revision is kept to freeze the source code to the time of manuscript submission.

## Dependencies

We recommend using the Singularity image provided on Singularity Hub or Sylabs Cloud,
which contains a conda environment with all CANDIA dependencies.
Alternatively, the image may be built from the supplied `candia.def` file 
(requires root permissions).

To download the image from either location:

* Singularity Hub: `singularity pull shub://fburic/candia:latest`
* Sylabs Cloud: `singularity pull library://fburic/default/candia` 

A conda environment can also be built from scratch using the provided `candia_env.yaml`
specification file. Here, commands will be shown using the Singularity image.

## Usage

The pipeline currently consists of a collection of scripts.
Step-by-step instruction are listed here using the toy data files provided in this repo. 
These are primarily meant to test the pipeline is working properly.

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

The execution of the pipeline is configured through a YAML file
(`test_experiment/config/candia.yaml`). 
This configuration file specifies the location of input and intermediate files,
as well as algorithm parameters. The paths are interpreted as relative to the 
experiment directory `root_dir`.

Please read through the next steps to see which parameters are relevant at each step.
A full configuraiton file is provided in the test experiment.

The scripts are expected to be run from the CANDIA top-level directory.

###  1. Convert DIA scan files from mzML to CSV

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

```bash
configfile='test/test_experiment/config/candia.yaml'

docker run --mount type=bind,source="$(pwd)"/test/test_experiment,target=/app/test/test_experiment candia:test \
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

```bash
docker run --mount type=bind,source="$(pwd)"/test/test_experiment,target=/app/test/test_experiment candia:test \
    snakemake --forceall -s scripts/util/adjust_swaths.Snakefile --configfile  ${configfile}
```


### 3. Split samples into slices

The input scans are partitioned into swaths and RT windows of width `window_size_sec`.
Note that this will create many small CSV files in `slices_location/swath=value/rt_window=value` subdirectories. 

> Expected running time: 30 min

Relevant pipeline config values:

* `samples_adjusted_swaths` - location of output adjusted CSV files
* `window_size_sec` - the width of the RT windows
* `slices_location` - location of output slices

```bash
docker run --mount type=bind,source="$(pwd)"/test/test_experiment,target=/app/test/test_experiment candia:test \
    python scripts/util/split_csv_maps_to_slices.py --config ${configfile}
```


## 4. Generate tensor files for all slices

The slices are converted into NumPy (swath, rt, sample) tensors stored as `.npy` files. 
The Snakefile will start the conversions for each slice in parallel and each such job may
be submitted independently as a HPC cluster job.

> Expected running time: 10 - 40 min (depending on available resources)

Relevant pipeline config values:

* `mass_tol_ppm` - The mass tolerance (in PPM) of the scan acquisition

The m/z precision is 10 decimals and 
m/z partitions with less than 5 time points in any sample are filtered out.

```bash
docker run --mount type=bind,source="$(pwd)"/test/test_experiment,target=/app/test/test_experiment candia:test \
    snakemake --jobs 4 -s scripts/util/generate_slice_tensors.Snakefile \
    -R generate_slice_tensor --forceall --rerun-incomplete --configfile ${configfile} 
```

On a slurm-managed cluster, this would be:

```bash
snakemake --jobs 200 -s scripts/util/generate_slice_tensors.Snakefile \
  -R generate_slice_tensor --forceall --rerun-incomplete --configfile ${configfile} \
  --cluster "sbatch -A {account.number} -t 06:00:00 --ntasks 6"
```


## 5. Run PARAFAC decomposition

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

### Workstation

The workstation command has the syntax:

`scripts/parafac/decompose_workstation.sh  EXPERIMENT_CONFIG_FILE N_PARALLEL_DECOMP_PER_GPU [SNAKEMAKE_OPTIONS]`

Example running with Singularity (with 2 parallel decompositions per GPU):

```bash
singularity exec --nv candia.sif /bin/bash -c "scripts/parafac/decompose_workstation.sh ${configfile} 2"
```

### HPC cluster

The cluster command has the syntax:

`scripts/parafac/decompose_cluster.sh ${snake_configfile} N_PARALLEL_DECOMP_PER_GPU`

And is designed to be run as an array job on slurm.
Example (with 6 parallel decompositions per GPU):

```bash
sbatch --array=0-1 --gres=gpu:1 --time 24:00:00 \
    scripts/parafac/decompose_cluster.sh "${configfile}" 6
```

If Singularity is not available, you can use the
`decomopse_cluster_no_singularity.sh` script. 


## 6. Index all PARAFAC models and components
 
This is a necessary step for downstream analysis. It assigns a unique ID to each
PARAFAC model and component, as well as record the filenames of all tensor models.
 
TODO: expand, clarify

The cluster command is;

```bash
sbatch --time 00:20:00 --ntasks 1 --cpus-per-task 30 \
    run_python.sh scripts/parafac/models.py -c ${configfile}
```


## 7. Select Best Models

TODO: expand, clarify

First, collect the time modes and measure unimodality fractions.

```bash
sbatch --ntasks 1 --cpus-per-task 32 \
    run_python.sh scripts/parafac/collect_time_mode_values.py -c ${configfile}
```

Then, select models:

```bash
scripts/parafac/select_best_models.sh -c ${configfile}
```


## 8. Identify proteins using Crux or MS-GF+

Configure which tool to use through the configuration file

TODO: expand, clarify

```bash
sbatch -t 04:00:00 \
    run_python.sh scripts/parafac/id_models_concat.py -c ${configfile}
```


## 9. Build library

TODO: expand, clarify

```bash
snakemake -s scripts/quantification/build_library.Snakefile --configfile ${configfile}
```


## 10. Quantify proteins with DIA-NN

TODO: expand, clarify

```bash
sbatch scripts/quantification/diann_slurm.sh --configfile ${configfile}
```


## 11. De novo sequencing with Novor and DeepNovo

Configure which tool to use through the configuration file

TODO: expand, clarify

```bash
snakemake -s scripts/denovo/sequence_best_models.Snakefile --configfile ${configfile}
```
