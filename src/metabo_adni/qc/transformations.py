import numpy as np
import pandas as pd
from typing import Union
from metabo_adni.data import load


def impute_metabolites(dat_dict: dict[str, pd.DataFrame],
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
    '''
    print('=== Imputing metabolites ===')
    total_points_imputed = 0
    total_mets_imputed = []

    for key in dat_dict:
        metabo_names = load._get_metabo_col_names(dat_dict[key],
                                                  key)
        if platform == 'p180':
            dat = dat_dict[key].loc[dat_dict[key].index < 99999, metabo_names]
        elif platform == 'nmr':
            dat = dat_dict[key].loc[:, metabo_names]
        else:
            dat = []
            raise Exception('No valid metabolomics platform')
        mets_to_impute = dat.columns[dat.isna().any()]
        data_points_impute = dat.isna().sum().sum()
        total_mets_imputed.extend(mets_to_impute)
        total_points_imputed = total_points_imputed + data_points_impute
        print('We will impute ' +
              str(len(mets_to_impute)) +
              ' metabolites and ' +
              str(data_points_impute) +
              ' data points in the ' +
              str(key) +
              ' cohort.')
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
    return dat_dict
