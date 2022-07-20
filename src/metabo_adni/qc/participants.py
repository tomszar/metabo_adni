import numpy as np
import pandas as pd
from typing import Union
from metabo_adni.data import load


def remove_missing(dat_dict: dict[str, pd.DataFrame],
                   cutoff: float) -> dict[str, pd.DataFrame]:
    '''
    Remove participants from dataframes due to missing data greater than
    cutoff.

    Parameters
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe to modify.
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
        dat = dat_dict[key].loc[dat_dict[key].index < 99999, metabo_names]
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


def consolidate_replicates(dat_dict: dict[str, pd.DataFrame]) -> \
        dict[str, pd.DataFrame]:
    '''
    Consolidate replicates by estimating the average across replicates.

    Parameters
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe to modify.

    Returns
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe with participants
        consolidated.
    '''
    print('=== Consolidating replicates ===')
    for key in dat_dict:
        metabo_names = load._get_metabo_col_names(dat_dict[key],
                                                  key)
        dat = dat_dict[key].loc[dat_dict[key].index < 99999, metabo_names]
        duplicated_ID = dat.index[
            dat.index.duplicated()].unique()
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
    '''
    print('=== Removing non-fasting participants ===')
    fasting_dat = load.read_fasting_file(fasting_file)
    fasting_participants = fasting_dat[fasting_dat == 1].index
    for key in dat_dict:
        keep_participants = dat_dict[key].index.isin(
                                    fasting_participants)
        _print_removed(keep_participants, key)
        dat_dict[key] = dat_dict[key].loc[keep_participants]
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
        print('We will remove the following ' +
              str(to_remove) +
              ' participants in the ' +
              cohort + ' cohort:')
        if isinstance(participant_table, pd.DataFrame):
            print(participant_table)
