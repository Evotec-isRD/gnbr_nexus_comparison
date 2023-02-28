"""This script loads in a table from the Nexus database called 'protein-level-averages' and process the data
for instance by adding all know sysnonyms of chemicals, as well as 'exploding' chemical-protein pairs
 (i.e. one per row rather than aggregated)"""

import pandas_utilities as pu
import pandas as pd
import os
import time
from data_dependencies import DataDependencies
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

def explode_uniprots(df: pd.DataFrame) -> pd.DataFrame:
    """
    The input df has a chemical with the target proteins (UniProt ID) all grouped together in a str separated by a ";"
    This explodes the UniProt column.
    :param df: nexus data
    :return: data exploded on the uniprot column
    """
    df['UNIPROT_ACCESSION_NUMBERS'] = df['UNIPROT_ACCESSION_NUMBERS'].str.split(';')
    df = df.explode('UNIPROT_ACCESSION_NUMBERS')
    return df

def make_chem_name_lower_case(df: pd.DataFrame) -> pd.DataFrame:
    """
    :param df: nexus
    :return: makes target name column lower case
    """
    logger.info('making target name column lower case')
    df['chem_name'] = df['chem_name'].str.lower()
    return df

############# merging with other nexus tables ####################
def merge_nx_tables_for_chemNames_to_proteins(df: pd.DataFrame) -> pd.DataFrame:
    """
    This loads the file from the Oracle command which joins structure_form_ID to structure_ID
    so that all structures in nexus have chemical synonym annotations
    :param df:
    :return: merged df with nexusIDs and all known synonyms of that compound
    """
    df_des = pd.read_csv(DataDependencies.NEXUS_SYNONYMS) # df_des for descriptions/annotations
    df_des.columns = ['STRUCTURE_ID','chem_name']
    logger.info('merging protein_level_averages with annotation descriptions')
    start = time.time()
    df = pd.merge(df, df_des, on='STRUCTURE_ID', how='outer')
    logger.info('time taken to merge protein_level_averages with annotation descriptions = ', time.time() - start)
    return df

def group_on_strID_uniID(df: pd.DataFrame)-> pd.DataFrame:
    """
    :param df:
    :returns: Groupby structureID and uniprotID, and aggregate all chemical names into a single column
    """
    logger.info('groupby strutureID and uniproID, takes 80 seconds on laptop')
    start = time.time()
    df = (df
          .groupby(['STRUCTURE_ID','UNIPROT_ACCESSION_NUMBERS']) # NOTE:I tried grouping with target name and molecule function too, but gave 1.5 M rather than 2.5M so must be nesting the uniprotID which have common names.'TARGET_NAME', 'MOLECULE_FUNCTION'
          .agg({
                    'chem_name': lambda x: list(set(list(x))),
                    'TARGET_NAME': lambda x: list(set(list(x))),
                    'MOLECULE_FUNCTION': lambda x: list(set(list(x))),
                    'gnbr_name':lambda x: list(set(list(x)))

                })
          .reset_index()
          )
    logger.info('groupby strutureID and uniproID took', time.time()-start)

    return df


############# mapping functions ##################
def map_fuzzy_matched_annotation_names_to_gnbr(df: pd.DataFrame) -> pd.DataFrame:
    """
    This creates a column in the input df called 'gnbr_name' which maps from a dictionary - where gnbr chemicals have
    been matched with a 0.95 fuzzy matched threshold - to a nexus chemical name. This dictionary has been
    made using the python package 'rapidfuzz' that compares two strings.
    So an example in Nexus is 'protoporphyrin IX’ while in gnbr this could be 'protoporphyrin 9’
    :param df:
    :return: This creates a column in the input df called 'gnbr_name'
    """
    logger.info('matching fuzzy similar words from nexus to gnbr')

    fuzzMatch_nexus_gnbrWithThemesDrpDup = pd.read_pickle(DataDependencies.FUZZMATCHED_CHEMNAME_NEXUS_GNBR)
    ##### extract 1st tuple (chemical name) from what rapid fuzz returns
    fuzzMatch_nexus_gnbrWithThemesDrpDup['1st_tuple'] = fuzzMatch_nexus_gnbrWithThemesDrpDup.iloc[:, 2].apply(
        lambda x: x if x is None else x[0])
    ##### make everything lowercase just to stop any confusion later on with the dictionary mapping
    fuzzMatch_nexus_gnbrWithThemesDrpDup = pd.concat(
        [fuzzMatch_nexus_gnbrWithThemesDrpDup[col].astype(str).str.lower() for col in
         fuzzMatch_nexus_gnbrWithThemesDrpDup.columns], axis=1)
    ##### drop duplicates to speed up the process
    fuzzMatch_nexus_gnbrWithThemesDrpDup = fuzzMatch_nexus_gnbrWithThemesDrpDup.drop_duplicates(
        subset='First entity name, raw string', keep='last')
    ##### zip the matches up from Nexus map back to the original gnbr dataset
    dict_compare_gnbrdropdup_nexus = dict(zip(fuzzMatch_nexus_gnbrWithThemesDrpDup.iloc[:, 1], #keeping this here to remind myself which way round to do it comparing gnbr to nexus
                                              fuzzMatch_nexus_gnbrWithThemesDrpDup.iloc[:, 3]))
    dict_compare_nexus_gnbrdropdup = dict(zip(fuzzMatch_nexus_gnbrWithThemesDrpDup.iloc[:, 3],
                                              fuzzMatch_nexus_gnbrWithThemesDrpDup.iloc[:, 1]))
    ##### Map matches in the dic back to gnbr data
    df['gnbr_name'] = df['chem_name'].str.lower().map(dict_compare_nexus_gnbrdropdup)
    return df

def main() -> pd.DataFrame:
    """
    This loads the nexus table called 'protein level averages' - which contains compounds their protein targets -
    or it make it by joining nexus protein targets to their annotations (name synonyms)
    :returns: a datafame with nexus data with proteins levels joined to annotations
    """
    nexus_preprocess_filepath = r'./data/df_nx_protein_level_averages_with_chemical_names.pkl'
    if os.path.exists(nexus_preprocess_filepath):
        logger.info(f'loading preprocessed nexus data from {nexus_preprocess_filepath}')
        df = pd.read_pickle(nexus_preprocess_filepath)
        return df
    else:
        logger.info(f'Starting nexus preprocessing pipeline')
        ##### this has 2,462,326 rows
        logger.info(f'loading nexus data from ./data/df_nexus_protein_level_averages.pkl')
        df = pd.read_pickle(DataDependencies.NEXUS_UNPROCESSED)

        df = (df
              ##### this must be done in order
              .pipe(merge_nx_tables_for_chemNames_to_proteins)  # this takes csv file which had annoatations and merges to protein averages
              .pipe(pu.remove_leading_trailing_white_space)  # takes few seconds
              .pipe(make_chem_name_lower_case)
              .pipe(map_fuzzy_matched_annotation_names_to_gnbr)  # This creates a column in the input df called 'gnbr_name' with a name in gnbr which matches that in nexus with 0.95 threshold
              .pipe(group_on_strID_uniID)
              .pipe(explode_uniprots)
              .pipe(pu.remove_leading_trailing_white_space)
              )
        ##### save the processed gnbr to pickle
        # df.to_parquet(output) # doesn't like some of the chemical names use pickle instead
        logger.info(f'saving processed nexus data to {nexus_preprocess_filepath}')
        df.to_pickle(nexus_preprocess_filepath)
        return df

if __name__ == '__main__':
    df = main()

