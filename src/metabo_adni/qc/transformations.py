import numpy as np
import pandas as pd
import scipy.stats as stats
from typing import Union
from metabo_adni.data import load


def imputation(dat_dict: dict[str, pd.DataFrame],
               platform: str,
               lod_directory: Union[str, None] = None) ->\
        dict[str, pd.DataFrame]:
    '''
    Impute metabolite concentration values by half the minimum observed value.
    If the LOD file is passed for the p180 platform, the imputation is
    half of the LOD value.

    Parameters
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe to modify.
    platform: str
        Metabolomics platform to process. Only done in the p180.
    lod_directory: Union[str, None]
        Path of the files containing the LOD information for the p180.

    Returns
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe with metabolites imputed.
    '''
    print('=== Imputing metabolites ===')
    total_points_imputed = 0
    total_mets_imputed = []

    for key in dat_dict:
        metabo_names = load._get_metabo_col_names(dat_dict[key],
                                                  key)
        indices = load._get_data_indices(dat_dict[key], platform)
        dat = dat_dict[key].loc[indices, metabo_names]
        mets_to_impute = dat.columns[dat.isna().any()]
        data_points_impute = dat.isna().sum().sum()
        total_mets_imputed.extend(mets_to_impute)
        total_points_imputed = total_points_imputed + data_points_impute
        print(f'{len(mets_to_impute)} metabolites and {data_points_impute} ' +
              f'data points will be imputed in the {key} cohort.')
        for j in mets_to_impute:
            indices = dat.loc[
                dat[j].isna()].index
            if platform == 'p180' and lod_directory is not None:
                lod_files = load.read_lod_files(lod_directory)
                dat = dat_dict[key].loc[dat_dict[key].index < 99999,
                                        list(metabo_names) +
                                        ['Plate.Bar.Code']]
                barcode = pd.Series(dat.loc[indices, 'Plate.Bar.Code'])
                vals = []
                for bar in barcode:
                    met_lod = lod_files[key].loc[:, 'Plate.Bar.Code'] == bar
                    vals.append(lod_files[key].loc[
                        met_lod, j])
                dat_dict[key].loc[indices, j] = np.mean(vals) * 0.5
            else:
                half_min = dat.loc[:, j].min() / 2
                dat_dict[key].loc[indices, j] = half_min
    print('')
    return dat_dict


def log2(dat_dict: dict[str, pd.DataFrame],
         platform: str) -> dict[str, pd.DataFrame]:
    '''
    Transform metabolite concentration values to log2 values.
    Add a constant of 1 before log transformation.

    Parameters
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe to modify.
    platform: str
        Metabolomics platform to process.

    Returns
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe with metabolites log2
        transformed.
    '''
    print('=== Log2 transformation ===')
    for key in dat_dict:
        metabo_names = load._get_metabo_col_names(dat_dict[key],
                                                  key)
        indices = load._get_data_indices(dat_dict[key], platform)
        dat = dat_dict[key].loc[indices, metabo_names]
        log2_dat = np.log2(dat)
        dat_dict[key].loc[indices, metabo_names] = log2_dat
    print('')
    return dat_dict


def zscore(dat_dict: dict[str, pd.DataFrame],
           platform: str) -> dict[str, pd.DataFrame]:
    '''
    Transform metabolite concentration values to zscore values.

    Parameters
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe to modify.
    platform: str
        Metabolomics platform to process.

    Returns
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe with metabolites zscore
        transformed.
    '''
    print('=== Zscore transformation ===')
    for key in dat_dict:
        metabo_names = load._get_metabo_col_names(dat_dict[key],
                                                  key)
        indices = load._get_data_indices(dat_dict[key], platform)
        dat = dat_dict[key].loc[indices, metabo_names]
        zscore_dat = dat.apply(stats.zscore,
                               nan_policy='omit')
        dat_dict[key].loc[indices, metabo_names] = zscore_dat
    print('')
    return dat_dict


def winsorize(dat_dict: dict[str, pd.DataFrame],
              platform: str) -> dict[str, pd.DataFrame]:
    '''
    Replace extreme values, 3 std, with cap values.

    Parameters
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe to modify.
    platform: str
        Metabolomics platform to process.

    Returns
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe with metabolites values
        winsorized.
    '''
    print('=== Winsorize metabolite values ===')
    for key in dat_dict:
        metabo_names = load._get_metabo_col_names(dat_dict[key],
                                                  key)
        indices = load._get_data_indices(dat_dict[key], platform)
        dat = dat_dict[key].loc[indices, metabo_names]
        three_std = np.std(dat) * 3
        total_replacements = 0
        for i, std in enumerate(three_std):
            col = dat.columns[i]
            high = dat.loc[:, col] > std
            low = dat.loc[:, col] < -std
            dat.loc[high, col] = std
            dat.loc[low, col] = -std
            total_replacements += sum(high) + sum(low)
        print(f'Replaced {total_replacements} values in {key} cohort.')
        dat_dict[key].loc[indices, metabo_names] = dat
    print('')
    return dat_dict
