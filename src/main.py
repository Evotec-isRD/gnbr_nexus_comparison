"""
This is the main script which runs all the other modules as a pipeline to allow for all the data processing and
production of figures for the paper
"""
import nexus_preprocess
import gnbr_preprocess
import merge_nexus_gnbr
import values_for_venn_figures
import plotting_prec_recall_ROC
import plot_upset_plots
import make_tables


def main():
    ##### processing scripts #####
    nexus = nexus_preprocess.main()
    gnbr = gnbr_preprocess.main()
    nx_gnbr_merged = merge_nexus_gnbr.main(gnbr, nexus) # note : this is an inner merge

    ##### plotting scripts #####
    values_for_venn_figures.main() # Figures 1 and 2
    plot_upset_plots.main(nx_gnbr_merged) # Figure 3, 5 and SI Figure 3
    plotting_prec_recall_ROC.main(nx_gnbr_merged)  # Figure 4
    make_tables.main(nx_gnbr_merged) # Tables 1-3

if __name__ == "__main__":
    main()


