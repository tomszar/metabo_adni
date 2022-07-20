import os
import glob
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
        Dict of dataframe names and dataframes.
    '''
    os.chdir(directory)
    dir_files = glob.glob('*')
    platform_files = []
    if platform == 'p180':
        platform_files = {'ADNI1-UPLC': '',
                          'ADNI1-FIA': '',
                          'ADNI2GO-UPLC': '',
                          'ADNI2GO-FIA': ''}
        file_names = ['ADMCDUKEP180UPLC_01_15_16.csv',
                      'ADMCDUKEP180FIA_01_15_16.csv',
                      'ADMCDUKEP180UPLCADNI2GO.csv',
                      'ADMCDUKEP180FIAADNI2GO.csv']
        na_values = ['< LOD',
                     'No Interception',
                     '>Highest CS']
        index_cols = ['RID']

    elif platform == 'nmr':
        platform_files = {'NMR': ''}
        file_names = ['ADNINIGHTINGALE2.csv']
        na_values = ['TAG']
        index_cols = ['RID', 'VISCODE2']

    else:
        raise Exception('The platform should be p180 or nmr')

    for i, f in enumerate(file_names):
        if f in dir_files:
            key = list(platform_files.keys())[i]
            dat = pd.read_csv(f, na_values=na_values).\
                set_index(index_cols)
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
                         axis=0,
                         inplace=True)
        fasting_dat.loc[i] = val
    fasting_dat = fasting_dat.sort_index()
    return fasting_dat


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
    platform: str
        Metabolomics platform to process.

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
