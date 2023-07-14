"""
This is the main script which runs all the other modules as a pipeline to allow for all the data processing and
production of figures for the paper
"""
import pandas as pd

import nexus_preprocess
import gnbr_preprocess
import merge_nexus_gnbr
import values_for_venn_figures
import plotting_prec_recall_ROC
import plot_upset_plots
import make_tables

# # # run in pycharm console
# runfile(
#     "C:/gits/gnbr_nexus_comparison/"
#     "src/main.py",
#     wdir="C:/gits/gnbr_nexus_comparison",
# )

CHEMBL_SUBSET = True
GNBR_PREPROCESSED = True
MAKE_CHEMBL_GNBR_MERGE = False
CLEAN_SAVE_CHEMBL_GNBR_MERGE = False
def main():
    ##### processing scripts #####
    if CHEMBL_SUBSET:
        nx_gnbr_merged = pd.read_pickle('data/chembl_gnbr_merged.pkl')
    else:
        nexus = nexus_preprocess.main()
        if GNBR_PREPROCESSED:
            gnbr = pd.read_parquet('data/df_gnbr_withThemesmapped_selectedCols.parquet')
        else:
            gnbr = gnbr_preprocess.main()

        nx_gnbr_merged = merge_nexus_gnbr.main(gnbr, nexus) # note : this is an inner merge

    ##### plotting scripts #####
    if not CHEMBL_SUBSET:
        values_for_venn_figures.main() # Figures 1 and 2
        make_tables.main(nx_gnbr_merged)  # Tables 1-3

    plot_upset_plots.main(nx_gnbr_merged) # Figure 3, 5 and SI Figure 3
    plotting_prec_recall_ROC.main(nx_gnbr_merged)  # Figure 4

    if MAKE_CHEMBL_GNBR_MERGE:
        # this gets rows with a chemblID
        nx_gnbr_merged['molecule_chembl_id'] = nx_gnbr_merged['chem_str'].astype(str).str.extract('(chembl\d+)')
        nx_gnbr_merged = nx_gnbr_merged.dropna(subset=['molecule_chembl_id'])

        # upper case in preparation for merge
        nx_gnbr_merged['molecule_chembl_id'] = nx_gnbr_merged['molecule_chembl_id'].str.upper()
        nx_gnbr_merged['uniprotID'] = nx_gnbr_merged['uniprotID'].str.upper()

        # get chembl mechasnim of action
        chembl_mechs = pd.read_csv('data/chembl_action-type_all.csv')

        # merge
        chembl_gnbr_merged = nx_gnbr_merged.merge(chembl_mechs,
                                                  on=['uniprotID', 'molecule_chembl_id'],
                                                  how='right')

        chembl_gnbr_merged = chembl_gnbr_merged.dropna()

        chembl_gnbr_merged.to_pickle('data/chembl_gnbr_merged.pkl')

    if CLEAN_SAVE_CHEMBL_GNBR_MERGE: # CLEAN-UP to save as .csv for Zenodo

        chembl_gnbr_merged = pd.read_pickle('data/chembl_gnbr_merged.pkl')
        chembl_gnbr_merged = chembl_gnbr_merged.drop_duplicates(subset=['chemical, raw string', 'gene, raw string',
                                                                'Dependency path', 'Sentence, tokenized'])
        chembl_gnbr_merged = chembl_gnbr_merged.drop(columns=['sum', 'fuzz_match', 'chem_ID', 'protein_ID',
                                                              'theme_value', 'Unnamed: 0', 'mechanism_of_action',
                                                              'target_chembl_id', 'relationship'
                                                              ])
        chembl_gnbr_merged = chembl_gnbr_merged.rename(columns={'chem_str': 'chemical synonyms',
                                                                'action_type': 'chembl_relationship',
                                                                'molecule_chembl_id': 'chemical_chembl_id'
                                                                        }
                                                               )
        chembl_gnbr_merged = chembl_gnbr_merged.loc[:, ['PubMed ID', 'chemical, raw string', 'gene, raw string',
                                                'Dependency path', 'Sentence, tokenized',
                                                '(A+) agonism', '(A-) antagonism', '(B) binding', '(E) affects prod.',
                                                '(E+) incr. expr/prod', '(E-) decr. expr/prod', '(K) metabol/pharmk',
                                                '(N) inhibits', '(O) transport', '(Z) enzyme', 'chemical synonyms',
                                                        'chemical_chembl_id',
                                                        'uniprotID',
                                                'protein_name', 'chembl_relationship']]

        chembl_gnbr_merged.to_csv('data/chembl_gnbr_merge_20k_sentences.csv')


if __name__ == "__main__":
    main()


