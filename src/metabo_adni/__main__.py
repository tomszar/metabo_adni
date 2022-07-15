import os
import argparse
from .data import load
from .qc import metabolites


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
    parser.add_argument('-c', '--cv',
                        type=float,
                        default=0.2,
                        help='Select the cutoff to remove metabolites\
                              based on CV values\
                              Default: 0.2.')
    parser.add_argument('-i', '--icc',
                        type=float,
                        default=0.65,
                        help='Select the cutoff to remove metabolites\
                              based on ICC values\
                              Default: 0.65.')
    args = parser.parse_args()
    print(args.platform)
    files = load.read_files(args.directory,
                            args.platform)
    files_m = metabolites.remove_missing(files,
                                         args.missing)
    files_m_cp = metabolites.cross_plate_correction(files_m,
                                                    args.platform)
    files_m_cp_cv = metabolites.remove_cv(files_m_cp,
                                          args.platform,
                                          args.cv)
    files_m_cp_cv_icc = metabolites.remove_icc(files_m_cp_cv,
                                               args.platform,
                                               args.icc)
    print('')
    print('=== Saving cleaned files ===')
    for key in files_m_cp_cv_icc:
        files_m_cp[key].to_csv(key + '.csv')
