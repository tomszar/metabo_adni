import pandas as pd
from metabo_adni.data import load


def remove_missing(dat_dict: dict[str, pd.DataFrame],
                   cutoff: float) -> dict[str, pd.DataFrame]:
    '''
    Remove participants from dataframes due to missing data greater than cutoff.

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
        dat = dat_dict[key][metabo_names] # REMOVE POOL SAMPLES
        total_cols = dat.shape[1]
        remove_part = dat.isna().sum(axis=1) / total_cols > cutoff

    print('')
    return dat_dict