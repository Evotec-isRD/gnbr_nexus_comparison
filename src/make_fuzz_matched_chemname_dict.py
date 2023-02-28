"""This is a stand-alone script which outputs a dictionary of 'fuzz-matched' chemical names between nexus and GNBR at a
trheshold of 0.95. This is saved in r'./data/df_chemGeneSortedWithTheme_dropDup_NexusFuzzMatch_unzip.pkl' which is used
in nexus_preprocess.py. The reason we have included it as a data dependency rather than generating the file  - is that
generating the file takes a long time (about an hour using 80 cores on a HPC).
"""

import pandas as pd
from pandarallel import pandarallel
from rapidfuzz import process, fuzz

# load in nexus annotations database (~1GB)
df_nexus_annotations = pd.read_parquet(r'./data/pddf_nx_annotation_proteinLevelAverages_explode.parquet')

# join group annotions by nexus structure ID
df_nexus_group = df_nexus_annotations.groupby('STRUCTURE_FORM_ID')['ANNOTATION','DESCRIPTION'].agg(list).reset_index()

# nexus annotations as a list (there are >3.5M), for comparison with gnbr chemical names as a list used later on.
nexus_annotationsSet = set([item.lower() for sublist in df_nexus_group['ANNOTATION'] for item in sublist])
nexus_annotationsSetLst = list(nexus_annotationsSet)

# load in gnbr chem-gene dataset
df_chemGeneSortedWithTheme = pd.read_csv('part-ii-dependency-paths-chemical-gene-sorted-with-themes.txt.gz', header=None, sep='\t')

# add column names for gnbr dataset
columns = ['PubMed ID',
 'Sentence number (0 = title)',
 'First entity name, formatted',
 'First entity name, location (characters from start of abstract)',
 'Second entity name, formatted',
 'Second entity name, location',
 'First entity name, raw string',
 'Second entity name, raw string',
 'First entity name, database ID(s)',
 'Second entity name, database ID(s)',
 'First entity type (Chemical, Gene, Disease)',
 'Second entity type (Chemical, Gene, Disease)',
 'Dependency path',
 'Sentence, tokenized']
df_chemGeneSortedWithTheme.columns = columns

# get chemical names from gnbr as a list. there are 73,275 unique names
gnbr_chemicalNamesSet = set([x for x in df_chemGeneSortedWithTheme.iloc[:, 6].str.lower()])

# initialise padarellel for multi-CPU usage
pandarallel.initialize()

# define a fuzzy matching function that compare two names and returns a score
def fuzzMatch(x):
    try:
        match = process.extractOne(x, nexus_annotationsSetLst, scorer=fuzz.ratio, score_cutoff=90)
        return match
    except:
        return None

# drop duplicates in the gnbr dataset chemical names
df_chemGeneSortedWithTheme_dropDup = df_chemGeneSortedWithTheme['First entity name, raw string']\
                                            .drop_duplicates(keep='last').reset_index()

# add a new column called 'nexus_fuzz_match' which returns a match between nexus and gnbr if the two str are >0.95 similar
df_chemGeneSortedWithTheme_dropDup['nexus_fuzz_match'] = df_chemGeneSortedWithTheme_dropDup['First entity name, raw string']\
                                                                    .astype(str).parallel_apply(fuzzMatch)
# save as a pkl
df_chemGeneSortedWithTheme_dropDup.to_pickle(r'./data/df_chemGeneSortedWithTheme_dropDup_NexusFuzzMatch_unzip.pkl')