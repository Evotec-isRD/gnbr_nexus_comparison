"""
Loads gnbr chem-gene data if the pipeline has already been run and the file created.
If not the pipe will join the themes with the rest of the data, which in the original GNBR paper comes
as two separate csv. The loaded data here has 1,728,361 rows
"""

import json
import pandas as pd
import pandas_utilities as pu
import os
import time
import logging
from data_dependencies import DataDependencies
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


class GetColumnNames:
    """This is a list of the column names that are headers in the GNBR data but are not included in the csv files
    that are downloaded, so are added manually here"""
    logger.info('getting column names')
    column_names = ['PubMed ID',
                    'Sentence number (0 = title)',
                    'First entity name, formatted',
                    'First entity name, location (characters from start of abstract)',
                    'Second entity name, formatted',
                    'Second entity name, location',
                    # 'First entity name, raw string',
                    'chemical, raw string',
                    # 'Second entity name, raw string',
                    'gene, raw string',
                    # 'First entity name, database ID(s)',
                    'chemical, database ID(s)',
                    # 'Second entity name, database ID(s)',
                    'gene, database ID(s)',
                    'First entity type (Chemical, Gene, Disease)',
                    'Second entity type (Chemical, Gene, Disease)',
                    'Dependency path',
                    'Sentence, tokenized',
                    ]
class GetSelectColumns:
    """These are the columns of interest in the GNBR data"""
    logger.info('selecting column names')
    select_columns = ['PubMed ID',
                      'chemical, raw string',
                      'gene, raw string',
                      'chemical, database ID(s)',
                      'gene, database ID(s)',
                      # 'First entity type (Chemical, Gene, Disease)', # have left as a note for completeness
                      # 'Second entity type (Chemical, Gene, Disease)', # have left as a note for completeness
                      'Dependency path',
                      'Sentence, tokenized',
                      ]

def map_entrezGeneID_to_uniprotID(df: pd.DataFrame) -> pd.DataFrame:
    """
    Maps Entrez gene Ids to UniProt IDs, so they can be matched to Nexus entities. Uses a csv mapping file
    made by using the UniProt mapping API.
    :param df:
    :return: df with new column called 'uniprotID'
    """
    logger.info('mapping entrez ID to Uniprot ID')
    df_geneMapping =  pd.read_csv(DataDependencies.GENE_MAPPING)
    ##### convert P_ENTREZGENEID to integer
    df_geneMapping['P_ENTREZGENEID'] = df_geneMapping['P_ENTREZGENEID'].astype('Int64').astype(str)
    ##### check what datatype is in gnbr for geneID. it is a string then
    #df_chemGeneSortedWithTheme['Second entity name, database ID(s)'].dtype
    dict_entrezID_to_uniprot = dict(zip(df_geneMapping['P_ENTREZGENEID'], df_geneMapping['ACC']))
    df['uniprotID'] = df['gene, database ID(s)'].map(dict_entrezID_to_uniprot)
    return df

def map_entities_to_themes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Joins the theme values to the entities from GNBR (part i joining to part ii).
    :param df: gnbr chem-gene part ii
    :return: gnbr themes mapped to part ii
    """
    dict_theme_value = r'./data/dict_theme_value.json'
    if not os.path.exists(dict_theme_value):
        # print('mapping themes values to entities dependency path (takes about 4 minutes)')
        logger.info('logging: mapping themes values to entities dependency path (takes about 4 minutes)')
        df1 = pd.read_csv(DataDependencies.GNBR_DISTRIBUTIONS, sep='\t')
        start = time.time()
        theme_value_dict = {x['path']: x.drop('path').to_dict() for i, x in df1.iterrows()}
        # print('time to map them values to the dependency path = ', time.time() - start)
        logger.info('time to map them values to the dependency path = ', time.time() - start)
        df['theme_value'] = df.loc[:, 'Dependency path'].str.lower().map(theme_value_dict)
        json.dump(theme_value_dict, open(dict_theme_value, 'w' ))
        # theme_value_dict1 = pd.Series(theme_value_dict) # option to save as parquet not used
        # theme_value_dict.to_parquet(output)
        return df
    else:
        logger.info('opening theme_value_dictionary JSON (could take a minute)')
        theme_value_dict = json.load(open(dict_theme_value))
        df['theme_value'] = df.loc[:, 'Dependency path'].str.lower().map(theme_value_dict)
        return df


def main() -> pd.DataFrame:
    """
    Loads gnbr chem-gene data if the pipeline has already been run and the file created.
    If not the pipe will join the themes with the rest of the data, which in the original GNBR paper comes
    as two separate csv. The loaded data here has 1,728,361 rows
    :returns a df preprocessed GNBR data
    """
    gnbr_preprocess_filepath = r'./data/df_gnbr_withThemesmapped_selectedCols.parquet'
    if os.path.exists(gnbr_preprocess_filepath):
        logger.info(f'loading preprocessed gnbr data from {gnbr_preprocess_filepath}')
        df = pd.read_parquet(gnbr_preprocess_filepath)
        return df
    else:
        logger.info('cannot find GNBR data. Running pipline to create it . . .')
        column_names = GetColumnNames.column_names
        select_columns = GetSelectColumns.select_columns

        logger.info('logging: loading gnbr part ii chem gene data')
        df = pd.read_csv(DataDependencies.GNBR_CHEMGENES, header=None, sep='\t')

        df = (df
              ##### Do this in order #####
              .pipe(pu.add_column_names, column_names=column_names)
              .pipe(pu.columns_to_select, select_columns=select_columns)
              .pipe(map_entities_to_themes)
              .pipe(map_entrezGeneID_to_uniprotID)
              )

        ##### save the processed gnbr to parquet
        logger.info(f'saving preprocessed file to {gnbr_preprocess_filepath}')
        df.to_parquet(gnbr_preprocess_filepath)

if __name__ == '__main__':
    df = main()

