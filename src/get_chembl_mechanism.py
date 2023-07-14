from chembl_webresource_client.new_client import new_client
import pandas as pd
import re
import time

# run in pycharm console
# runfile(
#     "C:/gits/gnbr_nexus_comparison/"
#     "src/get_chembl_mechanism.py",
#     wdir="C:/gits/gnbr_nexus_comparison",
# )

chembl_ids = pd.read_csv('data/chembl_ids.csv')
if not chembl_ids:
        df = pd.read_pickle('data/chembl_vs_gnbr.pkl')
        df = df.explode('chem_name')
        df['chembl_id'] = df['chem_name'].str.extract('(chembl)', flags= re.IGNORECASE)
        df['chembl_id'] = df['chem_name'].str.startswith('chembl', na=False)
        df = df[df['chembl_id']]
        chembl_ids = df['chem_name'].drop_duplicates()
        chembl_ids.to_csv('data/chembl_ids.csv', index=False)

chembl_ids = [x.upper() for x in chembl_ids['chem_name']]
chembl_ids = list(chembl_ids)

df = pd.read_csv('data/chembl_action-type_all.csv')
if not df:
        start = time.time()
        begin = 0
        end = begin+100000
        dfs = []
        while end < 900000:
                try:
                        test = chembl_ids[begin:end]+['CHEMBL25', 'CHEMBL109480', 'CHEMBL846', 'CHEMBL109480']
                        begin = end
                        end = begin+100000
                        mechanism = new_client.mechanism
                        mech = mechanism.filter(molecule_chembl_id__in=test)\
                                .only(['molecule_chembl_id', 'action_type', 'mechanism_of_action', 'target_chembl_id'])
                        df = pd.DataFrame(list(mech))
                        dfs.append(df)
                        # df.to_csv('data/chembl_action-type_all.csv')
                        stop = time.time() - start
                except:
                        test = chembl_ids[800000:-1]+['CHEMBL25', 'CHEMBL109480', 'CHEMBL846', 'CHEMBL109480']
                        mechanism = new_client.mechanism
                        mech = mechanism.filter(molecule_chembl_id__in=test)\
                                .only(['molecule_chembl_id', 'action_type', 'mechanism_of_action', 'target_chembl_id'])
                        df = pd.DataFrame(list(mech))
                        dfs.append(df)
                        # df.to_csv('data/chembl_action-type_all.csv')
                        stop = time.time() - start
        df = pd.concat(dfs)
        df.to_csv('data/chembl_action-type_all.csv')

df_targ = pd.read_csv('data/chembl_targets_with_uniprotID.csv')
if not df_targ:
        target_lst = df['target_chembl_id'].to_list()
        start = time.time()
        target = new_client.target
        targ = target.filter(target_chembl_id__in=target_lst)\
                .only(['target_chembl_id', 'target_components', 'pref_name'])
        df_targ = pd.DataFrame(list(targ))
        df_targ['uniprotID'] = df_targ['target_components'].apply(lambda x: x[0]['accession'] if bool(x) else None)
        stop = time.time() - start
        df_targ.to_csv('data/chembl_targets_with_uniprotID.csv')

target_chembl_unprot_dic = dict(zip(df_targ['target_chembl_id'], df_targ['uniprotID']))
df['uniprotID'] = df['target_chembl_id'].map(target_chembl_unprot_dic)
df.to_csv('data/chembl_action-type_all.csv')





