"""PLots all upset plots in the paper  - that is Figure 3 and 5 and suppl. info figure 2 and 3 """
import pandas as pd
from upsetplot import from_indicators
from upsetplot import UpSet
from upsetplot import plot
from matplotlib import pyplot as plt
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

def normalisation_fraction_across_all_themes_per_row(df: pd.DataFrame) -> pd.DataFrame:
    """
    :param df:
    :return: normalises values in the themes of GNBR across rows (i.e. each sentence is normalised to itself)
    """

    all_cols = ['(N) inhibits', '(A+) agonism', '(B) binding', '(A-) antagonism',
                       '(E+) incr. expr/prod', '(E-) decr. expr/prod', '(E) affects prod.',
                       '(O) transport', '(K) metabol/pharmk', '(Z) enzyme'
                       ]
    df["sum"] = df.loc[:, all_cols].sum(axis=1)
    df.loc[:, all_cols] = df.loc[:, all_cols].div(df["sum"], axis=0)
    print(f'normalisation:', df.loc[1:10, all_cols])
    return df

def get_aspirin_eg(df: pd.DataFrame) -> pd.DataFrame:
    """
    Gets data from the merged data relating to aspirin and cyclooxygenase
    :param df:
    :returns: sampled df with just aspirin and cyclooxygenase in it
    """
    asp_cyclo = df[df['chemical, raw string'].str.contains('aspirin', na=False, case=False) \
                        & df['gene, raw string'].str.contains('cyclooxygenase', na=False, case=False)
                        ]
    return asp_cyclo

def upset_plot_on_gnbr(df: pd.DataFrame, name: str =None, threshold: float =None) -> None:
    """
    Plots upsets plots with a threshold on gnbr data for a comparison with Nexus ground truth data
    :param df: input dataframe
    :param name: Category name (e.g. Antagonist)
    :param threshold: Threshold value to make the upset plot
    :return: an upset plot image
    """
    ##### RIGHT ORDER actually use these names as easier to read
    all_cols = [  # chemical-gene
        '(A+) agonism',
        '(A-) antagonism',
        '(B) binding',
        '(E) affects prod.',
        '(E+) incr. expr/prod',
        '(E-) decr. expr/prod',
        '(K) metabol/pharmk',
        # gene-chemical
        '(N) inhibits',
        '(O) transport',
        '(Z) enzyme']
    cols = all_cols
    df_up = df.loc[:, cols]
    ##### set a threshold
    threshold = threshold
    df_thr = df_up >= threshold
    ##### upset plot
    rel_types_gnbr = from_indicators(df_thr.columns, data=df_thr)
    upset = UpSet(rel_types_gnbr, show_percentages=True, show_counts=False)
    sub_cols = ['(N) inhibits', '(A+) agonism', '(B) binding', '(A-) antagonism']

    upset.style_subsets(present=sub_cols,
                        facecolor="blue",
                        ) #label="special"
    upset.plot()
    # plt.tight_layout()
    # plt.title('Relationships assigned to chemical-gene pairs with a \n threshold of {} on the EBC score'.format(threshold))
    plt.suptitle(
        'Relationships assigned to {} \n with a threshold of {} on the EBC score'.format(name, threshold))
    plt.show()
    plt.close()
    plt.savefig(f'./data/nexus_{name}_at_gnbr_threshold_{threshold}.png', format='png', dpi=300)

def main(nx_gnbr_merged: pd.DataFrame = None):
    """
    Runs scripts to plot upset plots
    :param nx_gnbr_merged: (this can either be passed in as a df, or it gets loaded from file if it already exists)
    :return: upset plots as matplotlib figures and save as png in data folder
    """
    if nx_gnbr_merged is None:
        logging.info('loading merged data')
        df = pd.read_pickle(r'./data/df_gnbr_with_themes_not_normalised_merged_to_nexus_169k.pkl')
    else:
        df = nx_gnbr_merged

    df = normalisation_fraction_across_all_themes_per_row(df)
    asp_cyclo = get_aspirin_eg(df)

    ##### plot upset section ######
    # on aspirin-cyclo example - this is Figure 3
    logging.info('plotting upset plots')
    upset_plot_on_gnbr(asp_cyclo, name='aspirin-cylcoxidase', threshold=0.00001) # change to 0.5 threshold for Figure 3

    # split data into the four Nexus classes, for Figure 5
    inh = df[df['relationship'].astype(str).str.contains('Inhibitor', na=False, case=False)]
    agon = df[df['relationship'].astype(str).str.contains('Agonist', na=False, case=False)]
    bind = df[df['relationship'].astype(str).str.contains('Binding', na=False, case=False)]
    ant = df[df['relationship'].astype(str).str.contains('Antagonist', na=False, case=False)]

    # plot upset for Figure 5
    upset_plot_on_gnbr(inh, name='Inhibitors', threshold=0.9)
    upset_plot_on_gnbr(agon, name='Agonists', threshold=0.9)
    upset_plot_on_gnbr(bind, name='Binding', threshold=0.9)
    upset_plot_on_gnbr(ant, name='Antagonists', threshold=0.9)

if __name__ == "__main__":
    main()