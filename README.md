# Metabo_ADNI

Metabolomics data processing for the ADNI data sets.
Currently, only supports the Biocrates p180 and Nightingale NMR platforms.

# Installation

- Clone the repo

```bash
git clone https://github.com/tomszar/adni_metabolomics.git
```

- Install metabo_adni

```bash
cd adni_metabolomics
pip install .
```

# Usage

In the folder with the required datasets, simply run:

```bash
clean_files
```

And metabo_adni will run with the default parameters.
**Note:** do not change the original name of the files.

## Options

- `-D`: define the directory were the files are located. Default, current working directory
- `-P`: define the platform, either p180 or nmr. Default, p180
- `-F`: define the fasting file. Default, BIOMARK.csv
- `--mmc`: remove metabolites with missing proportions greater than cutoff. Default, 0.2
- `--mpc`: remove participants with missing proportions greater than cutoff. Default, 0.2
- `--cv`: remove metabolites with CV values greater than cutoff. Default, 0.2
- `--icc`: remove metabolites with ICC values lower than cutoff. Default, 0.65
