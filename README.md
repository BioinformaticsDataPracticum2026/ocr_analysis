# ocr_analysis

Cross-species open chromatin region (OCR) comparison pipeline for the same organ across different organisms.

## Prerequisites

### Python Environment

```bash
module load python
python3 -m venv .venv
source .venv/bin/activate
```

### External Dependency

We use the `halLiftover-postprocessing` repository for `orthologFind.py`.

```bash
mkdir -p external
git clone https://github.com/pfenninglab/halLiftover-postprocessing.git external/halLiftover-postprocessing
```

## Run

```bash
python main.py
```

## Current Status

The project is still under development.

### Current priorities:

- keep the code simple
- make HALPER run correctly
- debug one step at a time before adding more pipeline stages

### Next Steps

- finish direct HALPER integration without relying on the shell wrapper
- generate orthologous OCR files for human-to-mouse and mouse-to-human
- add bedtools-based annotation for promoter/enhancer classification

## Contributions

Group Members: Qinglin Kong, Cindy Mai, Roy Hung, Shreya Nandakumar