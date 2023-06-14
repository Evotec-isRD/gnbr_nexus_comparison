import pandas as pd

df = pd.read_pickle('/Users/charliejeynes/PycharmProjects/gnbr_nexus_comparison/data/df_nx_protein_level_averages_with_chemical_names.pkl')
df = df.explode(column='chem_name')
df_chembl = df[df.chem_name.str.contains('chembl', na=False)]
df_chembl = df_chembl.explode(column=['TARGET_NAME'])
#
#
df_test = df.loc[0:100, :]
# df_test = df_test.explode(column='chem_name')
# df_chembl = df_test[df_test.chem_name.str.contains('chembl', na=False)]
# df_chembl_unq = df_chembl.drop_duplicates(subset=['chem_name', 'UNIPROT_ACCESSION_NUMBERS'])
df_chembl_unq = df_chembl.drop_duplicates(subset=['chem_name', 'TARGET_NAME'])
df_chembl_unq_sub = df_chembl_unq.loc[0:100, :]