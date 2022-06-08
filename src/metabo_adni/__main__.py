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
                              the ADNI files are located')
    parser.add_argument('-p', '--platform',
                        type=str,
                        default='p180',
                        help='Select the platform to analyze,\
                              either p180 or nmr')
    parser.add_argument('-m', '--missing',
                        type=float,
                        default=0.2,
                        help='Select the cutoff to remove metabolites\
                              based on missing data proportion')
    args = parser.parse_args()
    print(args.platform)
    files = load.read_files(args.directory,
                            args.platform)
    clean_files = qc.remove_missing_metabolites(files)
    for key in clean_files:
        clean_files[key].to_csv(key + '.csv')
