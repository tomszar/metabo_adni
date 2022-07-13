import pandas as pd
from metabo_adni.data import load


def remove_missing_metabolites(dat_dict: dict[str, pd.DataFrame],
                               cutoff: float) -> dict[str, pd.DataFrame]:
    '''
    Remove metabolites fram dataframes due to missing data greater than cutoff.

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
        Dictionary with dataframe name and dataframe with metabolites
        removed due to missingness.
    '''
    print('=== Removing metabolites with missing data greater than ' +
          str(cutoff) + ' ===')

    for key in dat_dict:
        col_names = load._get_metabo_col_names(dat_dict[key],
                                               key)
        dat = dat_dict[key][col_names]
        metabolite_table = _generate_missing_table(dat, cutoff)
        _print_metabolites_removed(metabolite_table, key)
        dat_dict[key].drop(metabolite_table.index,
                           axis=1,
                           inplace=True)
    return dat_dict


def remove_metabolites_cv(dat_dict: dict[str, pd.DataFrame],
                          platform: str,
                          cutoff: float) -> dict[str, pd.DataFrame]:
    '''
    Remove metabolites with a coefficient of variation (CV) higher than cutoff.
    The CV is computed taking into account duplicates and triplicates.
    Returns a cleaned dataframe.

    Parameters
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe to modify.
    platform: str
        Metabolomics platform to process. Only done in the p180.
    cutoff: float
        CV data removal cutoff.

    Returns
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe with corrected values.
    '''
    if platform == 'p180':
        print('=== Removing metabolites with CV values greater than ' +
              str(cutoff) + ' ===')
        for key in dat_dict:
            metabo_names = load._get_metabo_col_names(dat_dict[key],
                                                      key)
            dat = dat_dict[key].loc[dat_dict[key].index < 99999, metabo_names]
            cv_interplate = []
            duplicated_ID = dat.index[
                dat.index.duplicated()].unique()
            for id in duplicated_ID:
                duplicated_dat = dat.loc[id, :]
                cv = duplicated_dat.std() / duplicated_dat.mean()
                cv_interplate.append(cv)
            cv_interplate = pd.DataFrame(pd.DataFrame(cv_interplate).mean(),
                                         columns=['CV'])
            remove_met_table = cv_interplate[cv_interplate['CV'] > cutoff]
            _print_metabolites_removed(remove_met_table, key)
            dat_dict[key].drop(remove_met_table.index,
                               axis=1,
                               inplace=True)
        return dat_dict
    else:
        raise Exception('The platform should be p180 only')


def compute_cross_plate_correction(dat_dict: dict[str, pd.DataFrame],
                                   platform: str) -> dict[str, pd.DataFrame]:
    '''
    Computes the cross-plate correction for each metabolite
    by estimating the average metabolite value across sample pools
    and within plates to generate a per plate correction.
    Returns a corrected dataframe.

    Parameters
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe to modify.
    platform: str
        Metabolomics platform to process. Only done in the p180.

    Returns
    ----------
    dat_dict: dict[str, pd.DataFrame]
        Dictionary with dataframe name and dataframe with corrected values.
    '''
    if platform == 'p180':
        print('=== Applying a cross-plate correction ===')
        for key in dat_dict:
            metabo_names = load._get_metabo_col_names(dat_dict[key],
                                                      key)
            col_names = list(metabo_names) + ['Plate.Bar.Code']
            dat = dat_dict[key][col_names]
            dat_pool = dat.loc[999999, ]
            global_qc_average = dat_pool.loc[:, metabo_names].mean()
            plates_ID = dat_pool['Plate.Bar.Code'].unique()
            for j in range(len(plates_ID)):
                plate_name = plates_ID[j]
                # Estimate the correction
                q = dat_pool.loc[
                    dat_pool['Plate.Bar.Code'] == plate_name,
                    metabo_names].mean()
                correction = q / global_qc_average
                # Apply the correction
                correct_rows = (dat['Plate.Bar.Code'] == plate_name) &\
                    (dat.index < 99999)
                new_dat = dat.loc[correct_rows, metabo_names] / correction
                dat_dict[key].loc[correct_rows, metabo_names] = new_dat
        return dat_dict
    else:
        raise Exception('The platform should be p180 only')


def _generate_missing_table(dat: pd.DataFrame,
                            cutoff: float) -> pd.DataFrame:
    '''
    Produce a table with metabolites names and missing proportion less than
    cutoff.

    Parameters
    ----------
    dat: pd.DataFrame
        Dataframe with metabolites only.
    cutoff: float
        Missing data removal cutoff.

    Returns
    ----------
    metabolite_table: pd.DataFrame
        Table with metabolite names and missing proportion
    '''
    total_rows = len(dat)
    total_metabolites = pd.DataFrame(
        dat.isna().sum(axis=0) / total_rows,
        columns=['Missing percentage'])
    less_than_cutoff = dat.isna().sum(axis=0) / total_rows \
        > cutoff
    metabolite_table = total_metabolites[less_than_cutoff]

    return metabolite_table


def _print_metabolites_removed(metabolite_table: pd.DataFrame,
                               cohort: str) -> None:
    '''
    Print the list of metabolites that will be removed, and the value
    associated with the decision (either missing, CV, or ICC).

    Parameters
    ----------
    metabolite_table: pd.DataFrame
        Dataframe with metabolites names and either missing, CV or ICC
        values.
    cohort: str
        Cohort to append in print statement.

    Returns
    ----------
    None
    '''
    if len(metabolite_table) == 0:
        print('None of the metabolites will be dropped in the ' +
              cohort + ' cohort.')
    else:
        print('We will remove the following ' +
              str(len(metabolite_table)) +
              ' metabolites in the ' +
              cohort + ' cohort:')
        print(metabolite_table)
