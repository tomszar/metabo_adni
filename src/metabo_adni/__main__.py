import os
import argparse
from .data import load


def main():
    '''
    The main routine
    '''
    parser = argparse.ArgumentParser(
        description='Clean ADNI metabolomics datasets')
    parser.add_argument('-d', '--directory',
                        type=str,
                        nargs=1,
                        default=os.getcwd(),
                        help='Select the directory where\
                              the ADNI files are located')
    args = parser.parse_args()
    files = load.read_files(args.directory)
    print(files)
