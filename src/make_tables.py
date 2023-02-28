"""
This script makes table 2 and 3 in the paper
NOTES:
load gnbr all mentions no normalisation merged with nexus aggregated
drop_dups goes from 160k down to 10k
this is because in the merge, there are a few rows with different nexus relationships
"""

import pandas as pd
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

def normalisation_fraction_across_all_themes_per_row(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalises values in the themes of GNBR across rows (i.e. each sentence is normalised to itself)
    :param df:
    :return: all themes values in the df (see all_cols) are normalised across rows to between 0-1
    """
    all_cols = ['(N) inhibits', '(A+) agonism', '(B) binding', '(A-) antagonism',
                       '(E+) incr. expr/prod', '(E-) decr. expr/prod', '(E) affects prod.',
                       '(O) transport', '(K) metabol/pharmk', '(Z) enzyme'
                       ]
    df["sum"] = df.loc[:, all_cols].sum(axis=1)
    df.loc[:, all_cols] = df.loc[:, all_cols].div(df["sum"], axis=0)
    return df

def drop_dups_get_protein_names_col(df: pd.DataFrame) -> pd.DataFrame:
    """
    Drops duplicates on the 'chem raw str' and 'uniprot' columns
    :param df:
    :return: drops duplicates in the df and adds a protein name for the uniprot ID
    """
    df1 = df.drop_duplicates(subset=['chemical, raw string', 'uniprotID']).reset_index()
    ##### select columns #####
    cols = ['chemical, raw string', 'protein_name', '(N) inhibits', '(A+) agonism', '(B) binding', '(A-) antagonism',
            '(E+) incr. expr/prod', '(E-) decr. expr/prod', '(E) affects prod.', '(O) transport',
            '(K) metabol/pharmk', '(Z) enzyme', 'relationship']
    df1 = df1.loc[:, cols]
    df1= df1.assign(protein_name=lambda x: [x[0] for x in df1['protein_name']])
    return df1

def count_relationship_types_from_nexus(df1: pd.DataFrame) -> pd.DataFrame:
    """
    The relationship column  (from Nexus) is a list. For the example table we want to filter on examples where there
    is just one type between a compound and protein (e.g. 'inhibits') . So this function takes each row, take the first
    item from the list (to get a string), splits it, then counts it. This is for later processing where we want to
    filter for just one type.
    :param df1:
    :return: a new column called 'count' with how many relation type there are between a compound and protein
    """
    df1['relationship'] = df1['relationship'].apply(lambda x: x[0]) # list to string
    df1['relationship'] = df1['relationship'].apply(lambda x: x.split() if x is not None else 'None') # string to list for counts
    df1['count'] = df1['relationship'].apply(lambda x: len(x)) # count
    df1['relationship'] = df1['relationship'].apply(lambda x: x[0]) # list to string for later processing
    return df1

def sample_above_threshold(df1: pd.DataFrame) -> pd.DataFrame:
    """
    Sample above a gnbr relationship threshold and only on nexus samples that have 1 label
    :param df1:
    :return: a df that is a sample of the input df
    """
    themes = ['(N) inhibits', '(A+) agonism', '(B) binding', '(A-) antagonism',
                '(E+) incr. expr/prod', '(E-) decr. expr/prod', '(E) affects prod.', '(O) transport',
                '(K) metabol/pharmk', '(Z) enzyme']
    dic_dfs = []
    for theme in themes:
        dic_dfs.append(df1[(df1[theme] >= 0.5) & (df1[theme] <= 0.95) & (df1['count'] ==1) \
                          & (df1['relationship'].isin(['Antagonist', 'Agonist', 'Inhibitor', 'Binding']))] # (~df['relationship'].isin(['Substrate', 'Cofactor']))]
                        .sample(n=1, random_state=42) \
                        .round(2)\
                        .drop(columns='count'))
    df_tb = pd.concat(dic_dfs)
    df_tb = df_tb.rename(columns={'relationship': 'Nexus'})
    return df_tb

def save_score_values(df_tb: pd.DataFrame) -> pd.DataFrame:
    """
     saves examples with relationship normalised scores numbers - table 2 in paper
    :param df_tb: the output of sample_above_threshold()
    :return: saves an excel spreadsheet into 'data' folder
    """
    output = r'./data/df_tb_gnbr_nx_one_example_from_each_theme_theme_normalised.xlsx'
    df_tb = df_tb.sort_index()
    logger.info('saves examples with relationship normalised scores numbers - table 2 in paper to'
                './data/df_tb_gnbr_nx_one_example_from_each_theme_theme_normalised.xlsx')
    df_tb.to_excel(output)

def save_sentences(df: pd.DataFrame, df_tb: pd.DataFrame) -> pd.DataFrame:
    """
    get the sentences and the dependency paths from the same index as save_score_values()
    :param df: merged nexus/gnbr data
    :param df_tb: the output of sample_above_threshold()
    :return: saves an excel spreadsheet of into 'data' folder
    """
    df1 = df.drop_duplicates(subset=['chemical, raw string', 'uniprotID']).reset_index()
    sent_paths = df1.loc[df.index & df_tb.index, ['PubMed ID', 'Sentence, tokenized', 'Dependency path', 'uniprotID', 'chemical, raw string', 'gene, raw string']]
    output = r'./data/df_tb_gnbr_nx_one_example_from_each_theme_theme_normalised_sentences.xlsx'
    logger.info('saving sentence (table 3) ./data/df_tb_gnbr_nx_one_example_from_each_theme_theme_normalised_sentences.xlsx')
    sent_paths.to_excel(output)

def main(nx_gnbr_merged: pd.DataFrame = None) -> pd.DataFrame:
    """
    main function which pipes all the other functions
    :param nx_gnbr_merged:
    :return: saved excel from save_score_values() and save_sentences() which make table 2 and 3 in the paper
    """
    if nx_gnbr_merged is None:
        logger.info('loading ./data/df_gnbr_with_themes_not_normalised_merged_to_nexus_169k.pkl')
        df = pd.read_pickle(r'./data/df_gnbr_with_themes_not_normalised_merged_to_nexus_169k.pkl')
    else:
        logger.info('merging data')
        df = nx_gnbr_merged

    logger.info('making table')
    ##### make tables for paper using above function pipeline #####
    df = normalisation_fraction_across_all_themes_per_row(df)
    df1 = drop_dups_get_protein_names_col(df)
    df1 = count_relationship_types_from_nexus(df1)
    df_tb = sample_above_threshold(df1)
    save_score_values(df_tb)
    save_sentences(df, df_tb)

if __name__ == '__main__':
    main()