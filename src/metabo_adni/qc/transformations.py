import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
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


def merge(dat_dict: dict[str, pd.DataFrame],
          platform: str) -> dict[str, pd.DataFrame]:
    '''
    Merge dataframes across cohorts in the p180 platform.
    This action will delete all extra columns except for metabolites and
    index (RID), and all non-participant rows.
    In the nmr platform only deletes the columns.

    Parameters
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe to modify.
    platform: str
        Metabolomics platform to process.

    Returns
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Merged dataframe.
    '''
    print('=== Merging dataframes and leaving only metabolite columns and ' +
          'participants ===')
    dat_copy = dat_dict.copy()
    for key in dat_copy:
        metabo_names = load._get_metabo_col_names(dat_copy[key],
                                                  key)
        indices = load._get_data_indices(dat_copy[key], platform)
        dat = dat_copy[key].loc[indices, metabo_names]
        dat_copy[key] = dat
    if platform == 'p180':
        merge1 = dat_copy['ADNI1-FIA'].merge(dat_copy['ADNI1-UPLC'],
                                             how='inner',
                                             on='RID',
                                             suffixes=('_1FIA', '_1UPLC'))
        merge2 = dat_copy['ADNI2GO-FIA'].merge(dat_copy['ADNI2GO-UPLC'],
                                               how='inner',
                                               on='RID',
                                               suffixes=('_2FIA', '_2UPLC'))
        dat = pd.concat([merge1, merge2], join='inner')
        dat_dict = {'P180': dat}
    elif platform == 'nmr':
        dat_dict = dat_copy
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
        log2_dat = np.log2(dat + 1)
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


def residualize_metabolites(dat_dict: dict[str, pd.DataFrame],
                            platform: str) -> dict[str, pd.DataFrame]:
    '''
    Replace metabolite concentration values for the residuals on medication
    intake

    Parameters
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe to modify.
    platform: str
        Metabolomics platform to process.
    '''
    print('=== Replacing metabolite concentration values for residuals from ' +
          'medications ===')
    meds = load.read_meds_file()
    for key in dat_dict:
        metabo_names = load._get_metabo_col_names(dat_dict[key],
                                                  key)
        indices = load._get_data_indices(dat_dict[key], platform)
        dat = dat_dict[key].loc[indices, metabo_names]
        residuals = _get_residuals(dat, meds)
        dat_dict[key].loc[indices, metabo_names] = residuals
    return dat_dict


def _get_residuals(outcomes,
                   predictors):
    '''
    Get residuals for each column in the outcome dataframe.
    Both outcomes and predictors need an 'RID' column as index.

    Parameters
    ----------
    outcomes: pd.DataFrame
        dataframe with the outcome variables
    predictors: pd.DataFrame
        dataframe with the predictor variables
    '''
    new_dat = pd.merge(outcomes,
                       predictors,
                       left_on='RID',
                       right_on='RID')
    Y = new_dat.loc[:, outcomes.columns]
    residuals = Y.copy()
    X = new_dat.loc[:, predictors.columns]
    # Remove meds with only zeros
    keep_meds = X.mean() > 0
    X = X.loc[:, keep_meds]
    for y in Y:
        results = sm.OLS(exog=X,
                         endog=Y[y]).fit()
        med_names = list(X.columns)
        n_significants = sum(results.pvalues < 0.05)
        n_not_significants = sum(results.pvalues > 0.05)
        while n_not_significants > 0:
            drop_med = results.pvalues[results.pvalues > 0.05].\
                sort_values(ascending=False).\
                index[0]
            med_names.remove(drop_med)
            if not med_names:
                print(f'No significant medications in {y}')
                break
            else:
                results = sm.OLS(exog=X[med_names],
                                 endog=Y[y]).fit()
                n_significants = sum(results.pvalues < 0.05)
                n_not_significants = sum(results.pvalues > 0.05)
        if n_significants > 0:
            print(f'There are significant medications in {y}')
            print(med_names)
            print('')
            residuals[y] = results.resid
    return(residuals)
