# ADNI metabolomics QC

Metabolomics data processing for the ADNI data sets.
Currently, only supports the Biocrates p180 and Nightingale NMR platforms.

## Installation

metabo_adni is distributed as a python package, so install it by running:

```bash
pip install metabo_adni
```

## Usage

In the folder with the required datasets, simply run:

```bash
clean_files
```

And metabo_adni will run with the default parameters.
**Note:** do not change the original name of the files.

## Commands

- `-D`: define the directory were the files are located. Default, current working directory
- `-P`: define the platform, either p180 or nmr. Default, p180
- `-F`: define the fasting file. Default, BIOMARK.csv
- `-L`: define the directory were the LOD p180 files are located. Default, current working directory
- `--mmc`: remove metabolites with missing proportions greater than cutoff. Default, 0.2
- `--mpc`: remove participants with missing proportions greater than cutoff. Default, 0.2
- `--cv`: remove metabolites with CV values greater than cutoff. Default, 0.2
- `--icc`: remove metabolites with ICC values lower than cutoff. Default, 0.65
- `--log2`: apply log2 transformation to metabolite concentration values
- `--merge`: merge data frames across cohorts
- `--zscore`: apply zscore transformation to metabolite concentration values
- `--winsorize`: winsorize extreme values (more than 3 std of mean)
- `--remove-moutliers`: remove multivariate outliers using the Mahalanobis distance
- `--residualize-meds`: replace metabolite values with residuals from a regression with medication intake. Note that residuals are scaled to unit variance

## Files

The following are the files needed to run metabo_adni.

### NMR

For Nightingale's NMR metabolomics platform the file recommended is the one reanalyzed using the "2020 update", an updated quantification library.
The name of the item to download is `ADMC Nightingale Platform NMR Post-Unblinding Re-Analysis of Lipoproteins and Metabolites [ADNI1,GO,2]`, which downloads the `ADNINIGHTINGALE2.csv` file.

### Biocrates p180

For the Biocrates p180 platform, several files need to be downloaded that contain the proper data, QC tags, and LOD values:

#### Data

Four datasets contain the proper data, divided by method (FIA, UPLC) and cohort (ADNI1, ADNI2GO):

- `ADMCDUKEP180FIA_01_15_16.csv` obtained from `ADMC Duke Biocrates P180 Kit Flow injection analysis [ADNI1]` item
- `ADMCDUKEP180FIAADNI2GO.csv` obtained from `ADMC Duke Biocrates p180 Kit Flow injection analysis [ADNIGO,2]` item
- `ADMCDUKEP180UPLC_01_15_16.csv` obtained from `ADMC Duke Biocrates P180 Kit Ultra Performance Liquid Chromatography [ADNI1]` item
- `ADMCDUKEP180UPLCADNI2GO.csv` obtained from `ADMC Duke p180 Ultra Performance Liquid Chromatography [ADNIGO,2]` item

#### LOD

LOD values for the ADNI2GO can be found in the supplementary material `ADMC Duke p180 Supplementary Files [ADNIGO,2]`.
LOD values for the ADNI1 cohort coming soon.

### Fasting

To remove non-fasting participants, we need the `BIOMARK.csv` file downloaded from the `Biomarker Samples [ADNI1,GO,2,3,4]` item from ADNI.

### Medications

Medication intake information can be found in the `ADMC Duke ADNI2/GO Drug Classes` item, which donwloads the `ADMCPATIENTDRUGCLASSES_20170512.csv` file
