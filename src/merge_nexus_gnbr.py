"""
This script does processing of nexus and gnbr data to make them compatible for merging,
so that comparisons can be made.
This script has to be called from main.py in the src folder, as it's inputs are the preprocessed nexus and
gnbr data. It  does an inner merge on the two datasets on fuzzy matched chemical names and UniProt ID
"""
import pandas as pd
import pandas_utilities as pu
import os
import logging
################################ logging handling #####################################
log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logger = logging.getLogger(__name__)
# To override the default severity of logging
logger.setLevel('DEBUG')
# Use FileHandler() to log to a file
file_handler = logging.StreamHandler()
formatter = logging.Formatter(log_format)
file_handler.setFormatter(formatter)
# Don't forget to add the file handler
logger.addHandler(file_handler)

def gnbr_clean_theme_names(df_gnbr: pd.DataFrame) -> pd.DataFrame:
    """
    :param df_gnbr:
    :return: Renames the 'themes' columns in the gnbr dataset to something more human readable and also exclues the '.ind'
    columns
    """

    #####  these names as easier to read
    theme_short = [
        # chemical-gene
        '(A+) agonism',
        '(A-) antagonism',
        '(B) binding',
        '(E) affects prod.',
        '(E+) incr. expr/prod',
        '(E-) decr. expr/prod',
        '(K) metabol/pharmk',
        # gene-chemical
        '(N) inhibits',
        '(O) transport',
        '(Z) enzyme'
    ]

    ##### extract each theme from the dictionary in the column 'theme_value'
    row_themes = [x for x in df_gnbr['theme_value']]
    df_themes = pd.DataFrame(row_themes)
    get_cols = [col for col in df_themes.columns if '.ind' not in col]
    df_themes = df_themes.loc[:, get_cols]
    df_themes.columns = theme_short
    df_gnbr = df_gnbr.join(df_themes)

    return df_gnbr

def gnbr_remove_none_rows(df_gnbr: pd.DataFrame) -> pd.DataFrame:
    """
    :param df_gnbr:
    :return: remove None rows from gnbr to stop merge crashing in ID columns
     then remove None in uniprotIDs
    """
    df_gnbr = df_gnbr.drop(['chemical, database ID(s)', 'gene, database ID(s)'], axis=1).dropna()
    return df_gnbr

def nx_rename_columns(df_nx: pd.DataFrame) -> pd.DataFrame:
    """
    :param df_nx:
    :return: change name of columns add a new column as the df name
    """
    nx_col_names = ['chem_ID', 'protein_ID', 'chem_str', 'protein_name', 'relationship', 'fuzz_match']
    df_nx.columns = nx_col_names
    return df_nx

def nx_remove_synonymns_not_in_gnbr(df_nx: pd.DataFrame) -> pd.DataFrame:
    """
    :param df_nx:
    :return: take out instances of synonyms in nexus that is not in gnbr
    """
    df_nx = df_nx.explode('fuzz_match')
    df_nx = df_nx[df_nx['fuzz_match'].notna()]  # these are matches to gnbr 518231, but with all synonyms in gnbr
    return df_nx

def my_merge(df_gnbr: pd.DataFrame, df_nx: pd.DataFrame) -> pd.DataFrame:
    """
    :param df_gnbr: gnbr preprocessed data
    :param df_nx: nexus preprocessed data
    :return: merge the full gnbr dataset with unnormalised themes to nexus
    on fuzzy matched chemical names and UniProt ID
    """
    # default pandas merge is inner
    df_mg = df_gnbr.merge(df_nx,
                          left_on=['chemical, raw string', 'uniprotID'],
                          right_on=['fuzz_match', 'protein_ID', ]
                          )
    return df_mg


def main(gnbr: pd.DataFrame, nexus: pd.DataFrame) -> pd.DataFrame:
    """
    The pipeline here is to make the gnbr and nexus dataframes compatible for merging renaming columns.
    The file is then saved into a data dir as df_gnbr_with_themes_not_normalised_merged_to_nexus_169k.pkl
    the 169 or 161k refers to the rows.
    :param: gnbr is gnbr data
    :param: nexus is nexus data
    :returns: merged nexus and gnbr data
    """

    nx_gnbr_merged_filepath = './data/df_gnbr_with_themes_not_normalised_merged_to_nexus_169k.pkl'

    if (os.path.exists(nx_gnbr_merged_filepath)): # for the Venn diagrams I want to merge outer
        logger.info('loading merged nexus and gnbr data from ..data/df_gnbr_with_themes_not_normalised_merged_to_nexus_169k.pkl')
        df_mg = pd.read_pickle('./data/df_gnbr_with_themes_not_normalised_merged_to_nexus_169k.pkl')
        return df_mg
    else:
        logger.info('starting merging pipeline')
        gnbr = (
                gnbr
                .pipe(gnbr_clean_theme_names)
                .pipe(pu.make_all_lower_case)
                .pipe(pu.remove_leading_trailing_white_space)
                .pipe(gnbr_remove_none_rows)
                )

        nexus = (
                nexus
                .pipe(nx_rename_columns)
                .pipe(nx_remove_synonymns_not_in_gnbr)
                .pipe(pu.make_all_lower_case)
                .pipe(pu.remove_leading_trailing_white_space)
                )
        df_mg = my_merge(gnbr, nexus)
        logger.info(f'saving merged data to {nx_gnbr_merged_filepath}')
        df_mg.to_pickle(nx_gnbr_merged_filepath)
        return df_mg

if __name__ == "__main__":
    df_mg = main()