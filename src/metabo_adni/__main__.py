import os
import argparse
from .data import load
from .qc import metabolites
from .qc import participants


def main():
    '''
    The main routine
    '''
    parser = argparse.ArgumentParser(
        description='Clean ADNI metabolomics datasets')
    parser.add_argument('-D',
                        type=str,
                        default=os.getcwd(),
                        metavar='DIRECTORY',
                        help='Select the directory where\
                              the ADNI files are located.\
                              Default: current working directory.')
    parser.add_argument('-P',
                        type=str,
                        default='p180',
                        metavar='PLATFORM',
                        help='Select the platform to analyze,\
                              either p180 or nmr.\
                              Default: p180.')
    parser.add_argument('-F',
                        type=str,
                        default='BIOMARK.csv',
                        metavar='FASTING FILE',
                        help='Select the file from which to use the\
                              fasting information to keep only\
                              fasting participants. Default: BIOMARK.csv')
    parser.add_argument('--mmc',
                        type=float,
                        default=0.2,
                        metavar='MISSING METABOLITES CUTOFF',
                        help='Select the cutoff to remove metabolites\
                              based on missing data proportion.\
                              Default: 0.2.')
    parser.add_argument('--mpc',
                        type=float,
                        default=0.2,
                        metavar='MISSING PARTICIPANTS CUTOFF',
                        help='Select the cutoff to remove participants\
                              based on missing data proportion.\
                              Default: 0.2.')
    parser.add_argument('--cv',
                        type=float,
                        default=0.2,
                        metavar='COEFFICIENT OF VARIATION CUTOFF',
                        help='Select the cutoff to remove metabolites\
                              based on CV values.\
                              Default: 0.2.')
    parser.add_argument('--icc',
                        type=float,
                        default=0.65,
                        metavar='INTRA-CLASS CORRELATION CUTOFF',
                        help='Select the cutoff to remove metabolites\
                              based on ICC values.\
                              Default: 0.65.')
    args = parser.parse_args()
    files = load.read_files(args.D,
                            args.P)
    files = metabolites.remove_missing(files,
                                       args.mmc)
    files = metabolites.cross_plate_correction(files,
                                               args.P)
    files = metabolites.remove_cv(files,
                                  args.P,
                                  args.cv)
    files = metabolites.remove_icc(files,
                                   args.P,
                                   args.icc)
    files = participants.remove_missing(files,
                                        args.mpc)
    files = participants.consolidate_replicates(files)
    files = participants.remove_non_fasters(files,
                                            args.F)
    print('')
    print('=== Saving cleaned files ===')
    for key in files:
        files[key].to_csv(key + '.csv')
