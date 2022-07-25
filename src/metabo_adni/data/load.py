import os
import glob
import numpy as np
import pandas as pd
from typing import Union


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
        old_columns = dat.columns
        new_columns = old_columns.str.replace(pat='-|:|\(|\)| ',
                                              repl='.',
                                              regex=True)
        dat.columns = new_columns
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
    last_col = len(dat.columns)
    if 'ADNI1' in cohort:
        col_index = list(range(7, last_col-1))
    elif 'ADNI2GO' in cohort:
        col_index = list(range(24, last_col-1))
    elif 'NMR' in cohort:
        col_index = list(range(25, last_col-1))
    else:
        col_index = []

    col_names = dat.columns[col_index]
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


def check_p180_files(platform_files: dict[str, list[str]]) -> \
        dict[str, list[Union[str, bool]]]:
    '''
    Identify what types of file are there, and what type of QC can be done
    for the p180 platform

    Parameters
    ----------
    platform_files: dict[str, list[str]]
        Identified files in the directory

    Returns
    ----------
    p180_found_files: dict[str, list[str]]
        Types of QC available determined by file
    '''
    p180_needed_files = {'ADNI1-UPLC': ['ADMCDUKEP180UPLC_01_15_16.csv',
                                        'P180UPLCLODvalues_ADNI1.csv',
                                        '4610 UPLC p180 Data.xlsx'],
                         'ADNI1-FIA': ['ADMCDUKEP180FIA_01_15_16.csv',
                                       'P180FIALODvalues_ADNI1.csv',
                                       '4610 FIA p180 Data.xlsx'],
                         'ADNI2-UPLC': ['ADMCDUKEP180UPLCADNI2GO.csv',
                                        'P180UPLCLODvalues_ADNI2GO.csv',
                                        '4610 UPLC p180 Data.xlsx'],
                         'ADNI2-FIA': ['ADMCDUKEP180FIAADNI2GO.csv',
                                       'P180FIALODvalues_ADNI2GO.csv',
                                       '4610 FIA p180 Data.xlsx']}
    p180_found_files = {'ADNI1-UPLC': [False, False, False],
                        'ADNI1-FIA': [False, False, False],
                        'ADNI2-UPLC': [False, False, False],
                        'ADNI2-FIA': [False, False, False]}
    for k in p180_needed_files:
        for i, f in enumerate(p180_needed_files[k]):
            if f in platform_files['p180']:
                p180_found_files[k][i] = f
    return p180_found_files
