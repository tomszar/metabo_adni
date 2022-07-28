import os
import argparse
from .data import load
from .qc import metabolites
from .qc import participants
from .qc import transformations


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
    parser.add_argument('-L',
                        type=str,
                        default=os.getcwd(),
                        metavar='LOD FILES DIRECTORY',
                        help='Select the directory where\
                              the LOD p180 files are located.\
                              Default: current working directory.')
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
    parser.add_argument('--log2',
                        action='store_true',
                        help='Apply log2 transformation to metabolite\
                              concentration values.')
    parser.add_argument('--merge',
                        action='store_true',
                        help='Merge data frames across cohorts.\
                              Remove non-metabolite columns.')
    parser.add_argument('--zscore',
                        action='store_true',
                        help='Apply zscore transformation to metabolite\
                              concentration values.')
    parser.add_argument('--winsorize',
                        action='store_true',
                        help='Winsorize extreme values, i.e. more than 3 std,\
                              and replace them with 3 std')
    parser.add_argument('--remove-moutliers',
                        action='store_true',
                        help='Remove multivariate outliers using the\
                              Mahalanobis distance.')
    parser.add_argument('--residualize-meds',
                        action='store_true',
                        help='Replace metabolite concentration values\
                              with residuals from a regression on medication\
                              intake. Medication file needs to exist in the\
                              current directory.')
    args = parser.parse_args()
    files = load.read_files(args.D,
                            args.P)
    files = metabolites.remove_missing(files,
                                       args.P,
                                       args.mmc)
    if args.P == 'p180':
        files = metabolites.cross_plate_correction(files,
                                                   args.P)
        files = metabolites.remove_cv(files,
                                      args.P,
                                      args.cv)
        files = metabolites.remove_icc(files,
                                       args.P,
                                       args.icc)
    files = participants.remove_missing(files,
                                        args.P,
                                        args.mpc)
    files = participants.consolidate_replicates(files,
                                                args.P)
    files = participants.remove_non_fasters(files,
                                            args.F)
    files = participants.remove_bad_qc_tags(files,
                                            args.P)
    files = transformations.imputation(files,
                                       args.P,
                                       args.L)
    if args.log2:
        files = transformations.log2(files,
                                     args.P)
    if args.merge:
        files = transformations.merge(files,
                                      args.P)
    if args.zscore:
        files = transformations.zscore(files,
                                       args.P)
    if args.winsorize:
        files = transformations.winsorize(files,
                                          args.P)
    if args.remove_moutliers:
        files = participants.remove_moutliers(files,
                                              args.P)
    if args.residualize_meds:
        files = transformations.residualize_metabolites(files,
                                                        args.P)
    print('=== Saving cleaned files ===')
    for key in files:
        files[key].to_csv(key + '.csv')
