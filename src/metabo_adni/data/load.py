import os
import glob
from sys import platform
from typing import Union


def read_files(directory: str) -> dict[str, list[str]]:
    '''
    Read the files in the given directory and return a list of
    files with the possible parameters to run.

    Parameters
    ----------
    directory: str
        Directory where the files are located

    Returns
    ----------
    platform_files: dict[str, list[str]]
        Types of platforms found and their files associated
    '''
    os.chdir(directory)
    files = glob.glob('*')
    p180_files = []
    nmr_files = []

    for f in files:
        if 'ADMCDUKEP180' or 'p180 Data.xlsx' or 'LODvalues_ADNI' in f:
            p180_files.append(f)
        elif 'ADNINIGHTINGALE2' in f:
            nmr_files.append(f)

    platform_files = {'p180': p180_files,
                      'nmr': nmr_files}
    return(platform_files)


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
    return(p180_found_files)
