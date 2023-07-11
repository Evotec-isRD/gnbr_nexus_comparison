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

# run in pycharm console
# runfile(
#     "C:/gits/gnbr_nexus_comparison/"
#     "src/main.py",
#     wdir="C:/gits/gnbr_nexus_comparison",
# )


CHEMBL_SUBSET = True
GNBR_PREPROCESSED = True
GNBR_CHEMBL_INNER_MERGE = True
def main():
    ##### processing scripts #####

    if GNBR_CHEMBL_INNER_MERGE:
        nx_gnbr_merged = pd.read_pickle('data/df_gnbr_with_themes_not_normalised_merged_to_chembl_161k.pkl')
    else:
        if CHEMBL_SUBSET:
            nexus = pd.read_pickle('data/chembl_vs_gnbr.pkl')
        else:
            nexus = nexus_preprocess.main()

        if GNBR_PREPROCESSED:
            gnbr = pd.read_parquet('data/df_gnbr_withThemesmapped_selectedCols.parquet')
        else:
            gnbr = gnbr_preprocess.main()

        nx_gnbr_merged = merge_nexus_gnbr.main(gnbr, nexus, CHEMBL_SUBSET) # note : this is an inner merge


    ##### plotting scripts #####
    if not CHEMBL_SUBSET:
        values_for_venn_figures.main() # Figures 1 and 2

    plot_upset_plots.main(nx_gnbr_merged) # Figure 3, 5 and SI Figure 3
    plotting_prec_recall_ROC.main(nx_gnbr_merged)  # Figure 4
    make_tables.main(nx_gnbr_merged) # Tables 1-3

    if GNBR_CHEMBL_INNER_MERGE:
        nx_gnbr_merged = nx_gnbr_merged.drop(columns=['sum', 'fuzz_match', 'chem_ID', 'protein_ID'])
        nx_gnbr_merged = nx_gnbr_merged.rename(columns={'chem_str':'chemical synonymns',
                                                        'relationship': 'chembl_assigned_chem-prot_relationship'
                                                        }
                                                            )
        nx_gnbr_merged.to_csv('data/df_gnbr_with_themes_normalised_merged_to_chembl_161k.csv', index=False)

if __name__ == "__main__":
    main()


