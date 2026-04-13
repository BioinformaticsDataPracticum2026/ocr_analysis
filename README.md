# ocr_analysis

Cross-species open chromatin region (OCR) comparison pipeline for the same organ across different organisms.

## Prerequisites

This project is being developed on Bridges-2, a Linux HPC cluster.

### Recommended cluster setup

If this is your first time running the pipeline, create and activate a Python virtual environment:

```bash
module load python
module load bedtools
python3 -m venv .venv
source .venv/bin/activate
```

If you are in a fresh login session and the virtual environment already exists, run:

```bash
module load python
module load bedtools
source .venv/bin/activate
```

### External Dependencies

We use the `halLiftover-postprocessing` repository for `orthologFind.py`.

Clone it into `external/`:

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

Install the required packages, replacing `n` with the number returned by nproc:

```R
install.packages("BiocManager")
BiocManager::install(c("rGREAT", "GenomicRanges", "rtracklayer"), Ncpus = n)
```

When you are on a cluster, you may see a message saying that system library paths are not writable and asking whether you want to install packages into your home directory instead. Answer yes.


## Run

### Configurations

Before running the pipeline, you need to edit the configuration file. Configurations are located in `config.yaml`.

TODO: add how to configure `config.yaml`.

### Running the Pipeline

Running the pipeline is a simple one-liner, provided that your configurations are correct.

Because the first step (`halLiftover` / `HALPER`) can take a long time, the pipeline submits this task to the cluster queue. At the moment, this workflow is intended for SLURM-based clusters that support commands such as `squeue` and `scancel`.

```bash
python main.py
```

If you want to skip the first step and continue the pipeline using existing HALPER output, run:

```bash
python main.py --skip-halper
```

Even if you do not specify `--skip-halper`, the program can detect existing HALPER output and ask whether you want to reuse or overwrite those files.

If you need to check whether the job has been succesfully submitted, you can do so by using:

```bash
squeue -u $USER
```

## Current focus

This project is still under active development.

### Implemented

- `halLiftover/HALPER` to map mouse and human OCRs across species
- `bedtools intersect` to classify mapped OCRs as open or closed in the target species
- promoter/enhancer annotation of cross-species open OCRs using `bedtools closest` and target-species TSS annotations

### In Progres

- `rGREAT` analysis is working, but it still needs a Python wrapper

### Not Started

- motif analysis with `MEME` or `HOMER2`

## Contributions

QK (@umehina), CM (@cindymai-776), Roy Hung

## Special Acknowledgements

Thank you Isabella for pointing out how we should structure the project.