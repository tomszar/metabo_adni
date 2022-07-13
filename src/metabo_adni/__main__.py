import os
import argparse
import pandas as pd
from .data import load
from .qc import qc


def main():
    '''
    The main routine
    '''
    parser = argparse.ArgumentParser(
        description='Clean ADNI metabolomics datasets')
    parser.add_argument('-d', '--directory',
                        type=str,
                        default=os.getcwd(),
                        help='Select the directory where\
                              the ADNI files are located.\
                              Default: current working directory.')
    parser.add_argument('-p', '--platform',
                        type=str,
                        default='p180',
                        help='Select the platform to analyze,\
                              either p180 or nmr.\
                              Default: p180.')
    parser.add_argument('-m', '--missing',
                        type=float,
                        default=0.2,
                        help='Select the cutoff to remove metabolites\
                              based on missing data proportion\
                              Default: 0.2.')
    args = parser.parse_args()
    print(args.platform)
    files = load.read_files(args.directory,
                            args.platform)
    files_m = qc.remove_missing_metabolites(files)
    files_m_cp = qc.compute_cross_plate_correction(files_m,
                                                   args.platform)
    print('')
    print('=== Saving cleaned files ===')
    for key in files_m_cp:
        files_m_cp[key].to_csv(key + '.csv')
