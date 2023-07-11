import pandas as pd
import pickle

df = pd.read_pickle('data/df_nx_protein_level_averages_with_chemical_names.pkl')

def check_if_chembl_id(row: list) -> bool:
    bools = []
    for chem in row:
        chem = str(chem)
        if chem.find('chembl') != -1:
            bools.append(True)
        else:
            bools.append(False)
    return any(bools)

df['in_chembl'] = df.chem_name.apply(check_if_chembl_id)
df_chembl = df[df['in_chembl']]
df_chembl = df_chembl.drop(columns='in_chembl')
df_chembl.to_pickle('data/chembl_vs_gnbr.pkl')
