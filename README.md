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

- `-d`, `--directory`: define the directory were the files are located. Default, current working directory
- `-p`, `--platform`: define the platform, either p180 or nmr. Default, p180
- `-m`, `--missing`: remove metabolites with missing proportions greater than cutoff. Default, 0.2
