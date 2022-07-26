import warnings
import numpy as np
import pandas as pd
import scipy.stats as stats
from typing import Union
from metabo_adni.data import load


def remove_missing(dat_dict: dict[str, pd.DataFrame],
                   platform: str,
                   cutoff: float) -> dict[str, pd.DataFrame]:
    '''
    Remove participants from dataframes due to missing data greater than
    cutoff.

    Parameters
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe to modify.
    platform: str
        Metabolomics platform to process.
    cutoff: float
        Missing data removal cutoff.

    Returns
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe with participants
        removed due to missingness.
    '''
    print('=== Removing participants with missing data greater than ' +
          str(cutoff) + ' ===')

    for key in dat_dict:
        metabo_names = load._get_metabo_col_names(dat_dict[key],
                                                  key)
        indices = load._get_data_indices(dat_dict[key], platform)
        dat = dat_dict[key].loc[indices, metabo_names]
        total_cols = dat.shape[1]
        total_participants = pd.DataFrame(dat.isna().sum(axis=1) / total_cols,
                                          columns=['Missing percentage'])
        more_than_cutoff = dat.isna().sum(axis=1) / total_cols > cutoff
        participant_table = total_participants[more_than_cutoff]
        _print_removed(participant_table, key)
        if len(participant_table) > 0:
            dat_dict[key].drop(participant_table.index,
                               axis=0,
                               inplace=True)

    print('')
    return dat_dict


def consolidate_replicates(dat_dict: dict[str, pd.DataFrame],
                           platform: str) -> dict[str, pd.DataFrame]:
    '''
    Consolidate replicates by estimating the average across replicates.

    Parameters
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe to modify.
    platform: str
        Metabolomics platform to process.

    Returns
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe with participants
        consolidated.
    '''
    print('=== Consolidating replicates ===')
    # Ignore performance warnings from multiindex pandas
    warnings.simplefilter(action='ignore',
                          category=pd.errors.PerformanceWarning)
    for key in dat_dict:
        metabo_names = load._get_metabo_col_names(dat_dict[key],
                                                  key)
        indices = load._get_data_indices(dat_dict[key], platform)
        dat = dat_dict[key].loc[indices, metabo_names]
        duplicated_ID = dat.index[
            dat.index.duplicated()].unique()
        print(f'There are {len(duplicated_ID)} duplicated IDs in {key} cohort.')
        for j in duplicated_ID:
            # Averaging metabolites, and keeping first row of the other cols
            consolidated = dat.loc[j].mean()
            other_cols = dat_dict[key].columns.difference(metabo_names)
            first_row = dat_dict[key].loc[j, other_cols].iloc[0, ]
            new_row = pd.concat([first_row, consolidated])
            new_row = new_row[dat_dict[key].columns]
            dat_dict[key].drop(j,
                               axis='index',
                               inplace=True)
            dat_dict[key].loc[j] = new_row
        # The bl column is created in nmr (don't know why)
        if 'bl' in dat_dict[key].columns:
            dat_dict[key].drop('bl',
                               axis='columns',
                               inplace=True)
        dat_dict[key].sort_index(inplace=True)
    print('')
    return dat_dict


def remove_non_fasters(dat_dict: dict[str, pd.DataFrame],
                       fasting_file: str) -> dict[str, pd.DataFrame]:
    '''
    Remove participants that weren't fasting during data collection.

    Parameters
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe to modify.
    fasting_file: str
        Path of the file containing the fasting information.
        Minimum requiered columns are 'RID', and 'BIFAST'.

    Returns
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe with non-fasting
        participants removed.
    '''
    print('=== Removing non-fasting participants ===')
    fasting_dat = load.read_fasting_file(fasting_file)
    fasting_participants = fasting_dat[fasting_dat == 1].index
    for key in dat_dict:
        keep_participants = dat_dict[key].index.\
            get_level_values('RID').isin(
            fasting_participants)
        _print_removed(keep_participants, key)
        dat_dict[key] = dat_dict[key].loc[keep_participants]
    print('')
    return dat_dict


def remove_bad_qc_tags(dat_dict: dict[str, pd.DataFrame],
                       platform: str) ->\
        dict[str, pd.DataFrame]:
    '''
    Remove participants with bad QC tags.
    The QC tags in the nmr platform are in the main file.
    The QC tags in the p180 platform are in separate excel files,
    located in the qc_files_dir.

    Parameters
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe to modify.
    platform: str
        Metabolomics platform to process.

    Returns
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe with participants removed
        due to bad QC tags.
    '''
    print('=== Removing participants with bad QC tags ===')
    if platform == 'nmr':
        # Remove participants with at least one observed QC tag flagged.
        qc_tag_names = load._get_nmr_qc_cols()
        n_cols = len(qc_tag_names)
        id_list = dat_dict['NMR'].loc[:, qc_tag_names].iloc[:, :n_cols-1].\
            sum(axis=1) > 0
        to_remove = id_list.index[id_list].unique().get_level_values('RID')
        print('We will remove ' +
              str(len(to_remove)) +
              ' participants with bad QC tags in the nmr platform.')
        dat_dict['NMR'] = dat_dict['NMR'].drop(to_remove)
    elif platform == 'p180':
        # Only 'C5.DC..C6.OH.' in ADNI2GO FIA and
        # 'Taurine' in ADNI1 UPLC have bad QC flags
        for key in dat_dict:
            if key == 'ADNI1-UPLC':
                dat_dict[key].drop('Taurine',
                                   axis=1,
                                   inplace=True)
            elif key == 'ADNI2GO-FIA':
                dat_dict[key].drop('C5.DC..C6.OH.',
                                   axis=1,
                                   inplace=True)
    print('')
    return dat_dict


def remove_moutliers(dat_dict: dict[str, pd.DataFrame],
                     platform: str) -> dict[str, pd.DataFrame]:
    '''
    Remove multivariate outliers by computing the Mahalanobis distance

    Parameters
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe to modify.
    platform: str
        Metabolomics platform to process.

    Returns
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe with multivariate outliers
        removed.
    '''
    print('=== Removing multivariate outliers ===')
    for key in dat_dict:
        metabo_names = load._get_metabo_col_names(dat_dict[key],
                                                  key)
        indices = load._get_data_indices(dat_dict[key], platform)
        dat = dat_dict[key].loc[indices, metabo_names]
        cov_mat = dat.cov()
        p2 = dat.mean()
        cov_mat_pm1 = np.linalg.matrix_power(cov_mat, -1)
        distances = []
        for i, val in enumerate(dat.to_numpy()):
            p1 = val
            distance = (p1-p2).T.dot(cov_mat_pm1).dot(p1-p2)
            distances.append(distance)
        distances = np.array(distances)
        cutoff = stats.chi2.ppf(0.999, dat.shape[1]-1)
        i_to_remove = dat.loc[distances > cutoff].index
        print(f'{len(i_to_remove)} participants will be removed in the ' +
              f'{key} cohort')
        dat_dict[key].drop(i_to_remove,
                           axis='index',
                           inplace=True)

    print('')
    return dat_dict


def _print_removed(participant_table: Union[pd.DataFrame, np.ndarray],
                   cohort: str) -> None:
    '''
    Print the list of participants that will be removed, and the value
    associated with the decision (missing proportion).

    Parameters
    ----------
    participant_table: Union[pd.DataFrame, np.ndarray]
        Dataframe with participant names and missing proportion
        values (pd.DataFrame) or a list of which participants to keep
        (np.ndarray).
    cohort: str
        Cohort to append in print statement.

    Returns
    ----------
    None
    '''
    if len(participant_table) == 0:
        print('None of the participants will be dropped in the ' +
              cohort + ' cohort.')
    else:
        if isinstance(participant_table, pd.DataFrame):
            to_remove = len(participant_table)
        elif isinstance(participant_table, np.ndarray):
            to_remove = sum(~participant_table)
        else:
            to_remove = []
        print(f'{to_remove} participants will be removed in the {cohort} ' +
              'cohort.')
        if isinstance(participant_table, pd.DataFrame):
            print(participant_table)
