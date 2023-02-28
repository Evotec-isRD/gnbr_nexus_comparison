import pandas as pd
import numpy as np
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

def add_column_names(df: pd.DataFrame, column_names: list =None) -> pd.DataFrame:
    """
    Simple wrapper to add column names
    :param df:
    :param column_names:
    :return: df with column names added
    """
    if column_names is None:
        return df
    else:
        df.columns = column_names
        return df

def columns_to_select(df: pd.DataFrame, select_columns: list =None) -> pd.DataFrame:
    """
    Simple wrapper to select columns in a df
    :param df:
    :param select_columns:
    :return: smaller df with selected columns
    """
    logging.info(select_columns)
    if select_columns is None:
        return df
    else:
        df = df.loc[:, select_columns]
        return df

def make_all_lower_case(df: pd.DataFrame) -> pd.DataFrame:
    """
    Makes all str dtypes in a df lower case
    :param df:
    :return: lower cased strings
    """
    logging.info('making all strings lower case')
    # time comparison with applymap
    df = df.applymap(lambda s: s.lower() if type(s) == str else s)  # 17.7 seconds
    # with apply
    #     df = df.apply(lambda x: x.astype(str).str.lower()) # 19 seconds
    # with concat
    #     df  = pd.concat([df[col].astype(str).str.lower() for col in df.columns], axis=1) # 18.5 seconds
    return df

def remove_leading_trailing_white_space(df: pd.DataFrame) -> pd.DataFrame:
    """
    If the cell in a df is a str process it to remove leading and trailing white space from strings
    :param df:
    :return: str in df have leading and trailing white space from strings removed
    """
    print('removing leading and trailing white space from strings')
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    return df


def small_str_to_categorical(df: pd.DataFrame) -> pd.DataFrame:
    """
    Find str dtype in a dataframe and make categorical - this makes the df smaller in bytes
    :param df:
    :return: changes str to categorical dtype for smaller bytes
    """
    for name, column in df.head(100).iteritems():
        print(column.dtype)
        if column.dtype == object:
            print('yes')
            string_length_mean = np.mean([len(x) for x in column[column.notnull()]])
            print(string_length_mean)
            if string_length_mean < 20:
                # df = df.assign(name=lambda x: pd.Categorical(x[name]))
                df[name] = pd.Categorical(df[name])
    return df

# This is not used in any of the other scripts but is here for posterity
# def len_unique_values(df: pd.DataFrame, columns: list = None) -> pd.DataFrame:
#     """
#     Columns selects those that the user wants the number of unique values for.
#     If columns=None, it will automatically return the number of unique values in
#     Catergorical dtype columns.
#     :param df:
#     :param columns:
#     :return: a df with unique values from the input df
#     """
#     if columns is not None:
#         try:
#             df = df.iloc[:, columns]
#         except:
#             df = df.loc[:, columns]
#
#         num_of_unique = []
#         for name, col in df.iteritems():
#             num_of_unique.append({name: len(df[name].unique())})
#     else:
#         num_of_unique = []
#         for name, col in df.iteritems():
#             if col.dtype.name == 'category':
#                 num_of_unique.append({name: len(df[name].unique())})
#     dictionary = {k: v for dic in num_of_unique for k, v in dic.items()}
#     entities = {'entity': dictionary.keys(), 'number_uniques': dictionary.values()}
#     df_unq = pd.DataFrame(entities)
#     return df_unq

# BELOW FUNCTION IS A MAIN  BUT COMMENTING OUT BECAUSE WAS NEVER USED
# def readFilterTransform(filepath, header=None, sep='\t', column_names=None, select_columns=None) -> pd.DataFrame:
#     df = (pd.read_csv(filepath, header=None, sep='\t')
#           .pipe(add_column_names, column_names=None)
#           .pipe(columns_to_select, select_columns=None)
#           .pipe(make_all_lower_case)
#           .pipe(small_str_to_categorical)
#           )
#     return df





