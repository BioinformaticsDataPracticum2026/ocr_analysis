# ocr_analysis

Open chromatin region (OCR) comparison pipeline for the same organ across different organisms.

## Prerequisites

This project is being developed on Bridges-2, a Linux HPC cluster.

### Clone the repository
```bash
git clone https://github.com/BioinformaticsDataPracticum2026/ocr_analysis.git
```

### Recommended cluster setup

If this is your first time running the pipeline, create and activate a Python virtual environment _under the project directory_:

```bash
module load python bedtools homer
python3 -m venv .venv
source .venv/bin/activate
```

If you are in a fresh login session and the virtual environment already exists, run _under the project directory_:

```bash
module load python bedtools homer
source .venv/bin/activate
```

### halLiftover

We use `halLiftover`, which included in the HAL API. To install, follow the installation instuctions here.

https://github.com/ComparativeGenomicsToolkit/hal#installation

Before you click on that link, please understand that Linux clusters may already have HDF5 installed. To verify HDF5 installation, please use the following command

```bash
which h5cc
which h5dump
which h5ls
```

If one of those returns a path, HDF5 is likely available. If it is available, skip the section where it says `HDF5 1.10.1 with C++ API enabled`, and continue downloading and installing `sonLib`. You can skip the part where it says "optional support".

After building HAL, as mentioned in the HAL installation instructions, be sure to run the following command to add HAL to your `PATH` in the current shell:

```bash
export PATH=<path to hal>/bin:${PATH}
```

This makes it possible to use HAL from the terminal right away, without logging out and logging back in.

If you want this change to persist in future login sessions, add the same line to `~/.bashrc` by running:

```bash
echo 'export PATH=<path to hal>/bin:${PATH}' >> ~/.bashrc
source ~/.bashrc
```

Test with running one of the HAL API commands, such as `halLiftover`.

If you see the following `halLiftover` error, it means that HAL API has been successfully installed:
```bash
halLiftover
# Too few (required positional) arguments
# halLiftover v2.2: Map BED or PSL genome interval coordinates between two genomes.
# USAGE:
# halLiftover [Options] <halFile> <srcGenome> <srcBed> <tgtGenome> <tgtBed>
.....
```

### HALPER

We use the `halLiftover-postprocessing` repository for `orthologFind.py`.

When you are under the project root directory, clone it into `external/`:

```bash
mkdir -p external
git clone https://github.com/pfenninglab/halLiftover-postprocessing.git external/halLiftover-postprocessing
```

### R Dependency for rGREAT

For `rGREAT`, `GenomicRanges`, and `rtracklayer`, I recommend you to install packages in an interactive session on a node with enough memory and multiple CPU cores.

First, start an interactive session and check how many CPU cores are available:

```bash
interact
nproc
```

Then start R:

```bash
R
```

At least on Bridges-2, `R` is a pre-existing module, so you don't have to load it with `module load`. I also haven't figured out the module name for R, so thankfully it works out of the box.

In the R terminal, install the required packages, replacing `n` with the number returned by `nproc`:

```R
install.packages("BiocManager")
BiocManager::install(c("rGREAT", "GenomicRanges", "rtracklayer"), Ncpus = n)
```

When you are on a cluster, you may see a message saying that system library paths are not writable and asking whether you want to install packages into your home directory instead. Answer yes to those questions.

Quit R terminal with:

```R
q()
```

## Run

### Directory Structure

```
.
├── config.yaml
├── data -> /path/to/data
├── external
│   └── halLiftover-postprocessing
├── main.py
├── README.md
├── results
│   ├── bedtools
│   ├── halper
│   ├── rgreat
│   └── homer
├── scripts
│   ├── bedtools_preprocess.py
│   ├── bedtools_summary.py
│   ├── cross_species_ep.py
│   ├── halper.py
│   ├── homer.py
│   ├── open_closed.py
│   ├── promoter_enhancer.py
│   ├── __pycache__
│   ├── rgreat_online.R
│   └── rgreat.py
└── utils
    ├── check_dependencies.py
    ├── config.py
    ├── helpers.py
    └── __pycache__
```

- `config.yaml` the configuration file for the pipeline. The file contains info on how you would edit it.
- `external` contains external tools for this pipeline.
- `data` is a symlink to the actual location of the data. It is not neccessary to use symlinks, but we found it easier to be able to browse the raw data in VS Code while we were designing and implementing the pipeline.
- `main.py` is the entry point of the pipeline.
- `results` is the folder where all results are stored.
- `scripts` and `utils` contains source code to various parts of the pipeline.

**NOTE:** To set up symlinks to the data folder, you may run the following:

```bash
ln -s /path/to/data/source/ data
```

### Configurations

Before running the pipeline, you need to edit the configuration file. Configurations and instructions on how to edit these configuration are located in `config.yaml`.

### Running the Pipeline

Running the pipeline is a simple one-liner, provided that your configurations are correct. Run the following:

```bash
python main.py
```

Because the first step (`halLiftover` / `HALPER`) can take a long time, the pipeline submits this task to the cluster queue. At the moment, this workflow is intended for SLURM-based clusters that support commands such as `squeue` and `scancel`. `halLiftover` may take up to 24 hrs to run. Please refrain from running the pipeline again while you wait.

If you want to skip the first step and continue the pipeline using existing HALPER output, run:

```bash
python main.py --skip-halper
```

Even if you do not specify `--skip-halper`, the program can detect existing HALPER output and ask whether you want to reuse or overwrite those files.

You can also run rGREAT or Homer only by using the `--rgreat-only` or `--homer-only` flags.

```bash
python main.py --rgreat-only
python main.py --homer-only
```

### Managing Jobs Submitted to the Cluster

If you need to check whether the job has been succesfully submitted, you can do so by using:

```bash
squeue -u $USER
```

If you need to cancel all jobs, you can do so with:

```bash
scancel -u $USER
```

Please refer to the [SLURM documentation](https://slurm.schedmd.com/documentation.html) for additional commands. For PSC users, refer to the [user guide](https://www.psc.edu/resources/bridges-2/user-guide/).

## Contributions

QK (@umehina), CM (@cindymai-776), RH

## Citation

> **Kong, Q., Mai, C., & Hung, R. (2026).** *Open Chromatin Region Analysis for Human and Mice Pancreas Tissue.* 03-713: Bioinformatics Data Integration Practicum. Carnegie Mellon Univeristy. URL: https://github.com/BioinformaticsDataPracticum2026/ocr_analysis

## Special Acknowledgements

Thank you Isabella for pointing out how we should structure the project.