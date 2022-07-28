import os
import glob
import numpy as np
import pandas as pd


def read_files(directory: str,
               platform: str) -> dict[str, pd.DataFrame]:
    '''
    Read the files in the given directory and return a list of
    files with the possible parameters to run.

    Parameters
    ----------
    directory: str
        Directory where the files are located.
    platform: str
        Metabolomics platform to process.

    Returns
    ----------
    platform_files: dict[str, pd.DataFrame]
        Dictionary of dataframe names and dataframes.
    '''
    os.chdir(directory)
    dir_files = glob.glob('*')
    platform_files = []
    index_cols = ['RID']
    if platform == 'p180':
        df = pd.DataFrame()
        platform_files = {'ADNI1-UPLC': df,
                          'ADNI1-FIA': df,
                          'ADNI2GO-UPLC': df,
                          'ADNI2GO-FIA': df}
        file_names = ['ADMCDUKEP180UPLC_01_15_16.csv',
                      'ADMCDUKEP180FIA_01_15_16.csv',
                      'ADMCDUKEP180UPLCADNI2GO.csv',
                      'ADMCDUKEP180FIAADNI2GO.csv']
        na_values = ['< LOD',
                     'No Interception',
                     '>Highest CS']
    elif platform == 'nmr':
        platform_files = {'NMR': ''}
        file_names = ['ADNINIGHTINGALE2.csv']
        na_values = ['TAG']
    else:
        raise Exception('The platform should be p180 or nmr')

    for i, f in enumerate(file_names):
        if f in dir_files:
            key = list(platform_files.keys())[i]
            dat = pd.DataFrame(pd.read_csv(f,
                                           na_values=na_values)).\
                set_index(index_cols)
            dat = dat.sort_index()
            if 'ADNI2GO' in f:
                dat = _replace_bad_col_names(dat)
            # Carnosine is misspelled in ADNI2GO UPLC
            if f == 'ADMCDUKEP180UPLCADNI2GO.csv':
                dat = dat.rename(columns={'canosine': 'Carnosine'})
            platform_files[key] = dat

    return platform_files


def read_fasting_file(filepath: str) -> pd.Series:
    '''
    Read fasting file

    Parameters
    ----------
    filepath: str
        Path and name of the file.

    Returns
    ----------
    fasting_dat: pd.Series
        Series with fasting information.
    '''
    fasting_dat = pd.DataFrame(pd.read_csv(filepath,
                                           index_col='RID',
                                           na_values=-4))
    # Keep only information from baseline
    fasting_dat = fasting_dat.loc[fasting_dat.loc[:, 'VISCODE2'] == 'bl', ]
    fasting_dat = fasting_dat.loc[:, 'BIFAST']
    # If duplicates, keep the largest observed value
    duplicated_ID = fasting_dat.index[
        fasting_dat.index.duplicated()].unique()
    for i in duplicated_ID:
        val = fasting_dat.loc[i].max()
        fasting_dat.drop(i,
                         axis='index',
                         inplace=True)
        fasting_dat.loc[i] = val
    fasting_dat = fasting_dat.sort_index()
    return fasting_dat


def read_lod_files(directory: str) -> dict[str, pd.DataFrame]:
    '''
    Read the LOD files of the p180 platform.

    Parameters
    ----------
    directory: str
        Directory of the LOD files

    Returns
    ----------
    lod_files: dict[str, pd.DataFrame]
        Dictionary of LOD dataframe names and LOD dataframes.
    '''
    os.chdir(directory)
    lod_files = {'ADNI1-UPLC': '',
                 'ADNI1-FIA': '',
                 'ADNI2GO-UPLC': '',
                 'ADNI2GO-FIA': ''}
    filenames = ['P180UPLCLODvalues_ADNI1.csv',
                 'P180FIALODvalues_ADNI1.csv',
                 'P180UPLCLODvalues_ADNI2GO.csv',
                 'P180FIALODvalues_ADNI2GO.csv']
    for i, key in enumerate(lod_files):
        dat = pd.DataFrame(pd.read_csv(filenames[i],
                                       encoding='latin_1'))
        # Metabolite names in lod don't match those in the data
        # Replace '-', ':', '(', ')' and ' ' with '.'
        dat = _replace_bad_col_names(dat)
        if 'UPLC' in key:
            # Change metabolite name from Met.So to Met.SO
            dat.rename(columns={'Met.SO': 'Met.So'},
                       inplace=True)
        elif key == 'ADNI2GO-FIA':
            # In lod value ADNI2GO-FIA, the bar code plate needs fixing
            barcode = dat['Plate.Bar.Code']
            barcode = barcode.str.split(' ',
                                        expand=True)[2].\
                str.replace(pat='/',
                            repl='-')
            dat['Plate.Bar.Code'] = barcode
        lod_files[key] = dat
    return lod_files


def read_meds_file() -> pd.DataFrame:
    '''
    Reads the medication file.

    Returns
    ----------
    meds: pd.DataFrame
        Medication dataframe.
    '''
    file = 'ADMCPATIENTDRUGCLASSES_20170512.csv'
    file_exists = os.path.exists(file)
    if file_exists:
        meds = pd.read_csv(file).\
            set_index(['RID'])
        # Keeping only baseline
        baseline = meds.loc[:, 'VISCODE2'] == 'bl'
        # Removing extra column and no longer needed columns
        meds.drop(['NA', 'VISCODE2', 'Phase'],
                  axis='columns',
                  inplace=True)
        meds = meds.loc[baseline, :]
        # Replacing not NA values
        meds[meds.notna()] = 1
        # Replacing NA values
        meds = meds.replace(np.nan, 0)
    else:
        raise Exception('There is no medication file')
    return meds


def _get_metabo_col_names(dat: pd.DataFrame,
                          cohort: str) -> list[str]:
    '''
    Get the metabolite names based on position.

    Parameters
    ----------
    dat: pd.DataFrame
        Dataframe to retrieve column names.
    cohort: str
        The cohort of the dataframe (e.g. 'ADNI1-FIA').

    Returns
    ----------
    col_names: list[str]
        List of metabolite names
    '''
    cols = dat.columns
    if 'FIA' in cohort:
        start = list(cols).index('C0')
        end = list(cols).index('SM.C26.1') + 1
    elif 'UPLC' in cohort:
        start = list(cols).index('Ala')
        end = list(cols).index('SDMA') + 1
    elif 'NMR' in cohort:
        start = list(cols).index('TOTAL_C')
        end = list(cols).index('S_HDL_TG_PCT') + 1
    elif 'P180' in cohort:
        start = list(cols).index('C0')
        end = list(cols).index('SDMA') + 1
    else:
        start = 0
        end = 0

    col_names = cols[start:end]
    return col_names


def _get_data_indices(dat: pd.DataFrame,
                      platform: str) -> np.ndarray:
    '''
    Get the corresponding indices depending on the platform.

    Parameters
    ----------
    dat: pd.DataFrame
        Dataframe to retrieve column names.
    platform: str
        The platform of the dataframe ('p180' or 'nmr').

    Returns
    ----------
    indices: np.ndarray
        List of indices.
    '''
    if platform == 'p180':
        indices = dat.index < 99999
    elif platform == 'nmr':
        indices = np.repeat(True, len(dat.index))
    else:
        indices = np.ndarray([])
    return indices


def _get_nmr_qc_cols() -> list[str]:
    '''
    Get the QC column names for the nmr platforms

    Returns
    ----------
    qc_tag_names: list[str]
    '''
    qc_tag_names = ['EDTA_PLASMA',
                    'CITRATE_PLASMA',
                    'LOW_ETHANOL',
                    'MEDIUM_ETHANOL',
                    'HIGH_ETHANOL',
                    'ISOPROPYL_ALCOHOL',
                    'N_METHYL_2_PYRROLIDONE',
                    'POLYSACCHARIDES',
                    'AMINOCAPROIC_ACID',
                    'LOW_GLUCOSE',
                    'HIGH_LACTATE',
                    'HIGH_PYRUVATE',
                    'LOW_GLUTAMINE_OR_HIGH_GLUTAMATE',
                    'GLUCONOLACTONE',
                    'LOW_PROTEIN',
                    'UNEXPECTED_AMINO_ACID_SIGNALS',
                    'UNIDENTIFIED_MACROMOLECULES',
                    'UNIDENTIFIED_SMALL_MOLECULE_A',
                    'UNIDENTIFIED_SMALL_MOLECULE_B',
                    'UNIDENTIFIED_SMALL_MOLECULE_C',
                    'BELOW_LIMIT_OF_QUANTIFICATION']
    return qc_tag_names


def _replace_bad_col_names(dat: pd.DataFrame) -> pd.DataFrame:
    '''
    Replace columns with non-compatible names.

    Parameters
    ----------
    dat: pd.DataFrame
        Original dataframe with bad named columns.

    Returns
    ----------
    dat: pd.DataFrame
        Dataframe with replaced column names.
    '''
    old_columns = dat.columns
    new_columns = old_columns.str.replace(pat='-|:|\(|\)| ',
                                          repl='.',
                                          regex=True)
    dat.columns = new_columns
    return(dat)
