import pandas as pd
from metabo_adni.data import load


def remove_missing_metabolites(dat_dict: dict[str, pd.DataFrame],
                               cutoff: float = 0.2) -> None:
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
        print('None of the metabolites will be dropped')
    else:
        print('We will remove the following ' +
              str(len(metabolite_table)) +
              ' metabolites in the ' +
              cohort + ' cohort:')
        print(metabolite_table)
