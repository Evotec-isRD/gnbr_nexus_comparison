import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import PrecisionRecallDisplay
import seaborn as sns
import matplotlib.pyplot as plt
# sns.set(style='ticks', context='talk')
sns.set_theme(style="whitegrid")

import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
from sklearn.metrics import roc_curve, auc

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

def explode_drop_duplicates(df: pd.DataFrame) -> pd.DataFrame:
    """
    there are some duplications here in the gnbr dataset so this drops rows from 169073 to (this puzzles 111397)
    157646 sentence instances
    :param df:
    :return: df exploded on 'relationship' column
    """
    df1 = df.explode('relationship')
    df1 = df1.drop_duplicates(subset=['PubMed ID', 'chemical, raw string', 'gene, raw string',
                                        'Dependency path', 'Sentence, tokenized', 'uniprotID',
                                     'chem_ID', 'protein_ID', 'relationship'],
                             keep='last').reset_index()
    logger.info(f'experiment_drop_duplicates makes it goes from {len(df)} to {len(df1)}')
    return df1
def normalisation_fraction_across_all_themes_per_row(df: pd.DataFrame) -> pd.DataFrame:
    """
    normalises values in the themes of GNBR across rows (i.e. each sentence is normalised to itself)
    :param df:
    :return: noramlises values in the themes columns from 0-1
    """
    all_cols = ['(N) inhibits', '(A+) agonism', '(B) binding', '(A-) antagonism',
                       '(E+) incr. expr/prod', '(E-) decr. expr/prod', '(E) affects prod.',
                       '(O) transport', '(K) metabol/pharmk', '(Z) enzyme'
                       ]
    df["sum"] = df.loc[:, all_cols].sum(axis=1)
    df.loc[:, all_cols] = df.loc[:, all_cols].div(df["sum"], axis=0)
    return df
def mk_nx_ground_truth_labels(df: pd.DataFrame) -> pd.DataFrame:
    """
    This adds columns that normalise the relationship names between nexus and gnbr
    :param df:
    :return: harmonised relationship names between nexus and gnbr
    """
    df = df.loc[:, 'relationship'].reset_index().drop(columns='index')\
            .rename(columns={'relationship':'relationship_nexus'})
    ##### make sure it ignores 'inverse angonists/antagonist'
    df[df['relationship_nexus'].astype(str).str.contains('Inverse') == True] = 'ignore this cell as contained inverse'
    ##### contains
    df['(N) inhibits'] = df['relationship_nexus'].astype(str).str.contains('Inhibitor')
    df['(A+) agonism'] = df['relationship_nexus'].astype(str).str.contains('Agonist')
    df['(A-) antagonism'] = df['relationship_nexus'].astype(str).str.contains('Antagonist')
    df['(B) binding'] = df['relationship_nexus'].astype(str).str.contains('Binding')
    return df
def join_themes_gnbr_to_nexus(df: pd.DataFrame, df_nx: pd.DataFrame) -> pd.DataFrame:
    """
    This joins the theme scores from gnbr, to the binary classes in nexus ~160k rows
    :param df:
    :param df_nx:
    :return: df for comparison between nexus and gnbr
    """
    columns = ['(N) inhibits','(A+) agonism','(B) binding', '(A-) antagonism']
    df_gnbr = df.loc[:, columns]
    df_nx = df_nx.loc[:, columns]
    df_jn = df_gnbr.join(df_nx, lsuffix="_gnbr", rsuffix="_nx")
    return df_jn
def plot_counts_per_class_nx(df: pd.DataFrame) -> pd.DataFrame:
    """
    :param df:
    :return: Output is a bar chart which shows the total number of positive columns per class (inhibitor, agonist, binding,
    antagonist) and the total number of rows. This shows how balanced (or not) the data is.
    """
    cols = ['(N) inhibits_nx', '(A+) agonism_nx', '(B) binding_nx', '(A-) antagonism_nx']
    # n_cols = ['(N) inhibits', '(A+) agonism', '(B) binding', '(A-) antagonism']
    df1 = df.loc[:, cols]
    #
    n_cols = ['Inhibitor', 'Agonist', 'Binding', 'Antagonist', 'Total rows']
    df1['total'] = True
    # n_cols = ['Inhibitor', 'Agonist', 'Binding', 'Antagonist']
    df1.columns = n_cols # gets rid of the nx suffix
    df1 = pd.melt(df1, var_name='class', value_name='count')
    ##### add total observations bar
    # df2 = pd.DataFrame('class':np.pealen(df1))
    # df1 = df1.groupby('class')['count'].value_counts()
    df1 = df1[df1['count'] == True] # doesn't do anything
    print(f'true count in nexus is {len(df1)}') # here it selects all the true - which is the same as the at_least_one_true_in_nx
    ax = sns.countplot(y="class", data=df1, palette=["navy", "turquoise", "darkorange", "cornflowerblue", 'black'])
    plt.title('total number of positive labels in Nexus per class & total rows')
    plt.tight_layout()
    plt.show()
    plt.close()
    plt.savefig(f'./data/plots_counts_per_class_in_nexus.png', format='png', dpi=300)
    return ax
def get_yscore_ytest(df: pd.DataFrame) -> pd.DataFrame:
    """
    Preprocesses the data so that it is compatible with sci-kit learn
    :param df:
    :return: The output are two numpy arrays Y_test (gnbr data) to y_score (nexus labels) for comparison
    """
    nx_cols = ['(N) inhibits_nx', '(A+) agonism_nx', '(B) binding_nx', '(A-) antagonism_nx']
    Y_test = df.loc[:, nx_cols].to_numpy().astype(int)
    gnbr_cols = ['(N) inhibits_gnbr', '(A+) agonism_gnbr', '(B) binding_gnbr','(A-) antagonism_gnbr']
    y_score = df.loc[:, gnbr_cols].to_numpy()
    return Y_test, y_score

def plot_precision_recall(Y_test: np.array, y_score: np.array, normalisation:str ='per sentence') -> None:
    """
    :param Y_test:
    :param y_score:
    :param normalisation: # legacy from previous code where normalisation method was varied
    :return: plots and saves precision-recall curves
    """
    ##### get_micro_precision_recall():
    n_classes = Y_test.shape[1]
    # For each class
    precision = dict()
    recall = dict()
    average_precision = dict()
    for i in range(n_classes):
        precision[i], recall[i], _ = precision_recall_curve(Y_test[:, i], y_score[:, i])
        average_precision[i] = average_precision_score(Y_test[:, i], y_score[:, i])

    # A "micro-average": quantifying score on all classes jointly
    precision["micro"], recall["micro"], _ = precision_recall_curve(
                                                            Y_test.ravel(), y_score.ravel()
                                                        )
    average_precision["micro"] = average_precision_score(Y_test, y_score, average="micro")

    ##### plot_micro_averaged_precision_recall_curve():
    display = PrecisionRecallDisplay(
        recall=recall["micro"],
        precision=precision["micro"],
        average_precision=average_precision["micro"],
    )
    display.plot()
    _ = display.ax_.set_title("Micro-averaged over all classes")
    plt.show()

    ##### def plot_precision_recall_curve_for_each_class_and_isof1_curve():
    ##### setup plot details
    colors = cycle(["navy", "turquoise", "darkorange", "cornflowerblue", "teal"])
    _, ax = plt.subplots(figsize=(7, 8))
    f_scores = np.linspace(0.2, 0.8, num=4)
    lines, labels = [], []
    for f_score in f_scores:
        x = np.linspace(0.01, 1)
        y = f_score * x / (2 * x - f_score)
        (l,) = plt.plot(x[y >= 0], y[y >= 0], color="gray", alpha=0.2)
        plt.annotate("f1={0:0.1f}".format(f_score), xy=(0.9, y[45] + 0.02))
    display = PrecisionRecallDisplay(
        recall=recall["micro"],
        precision=precision["micro"],
        average_precision=average_precision["micro"],
    )
    display.plot(ax=ax, name="Micro-average precision-recall", color="gold")
    ##### I added classes name here
    classes = ['(N) inhibits', '(A+) agonism', '(B) binding', '(A-) antagonism']
    for i, color in zip(range(n_classes), colors):
        display = PrecisionRecallDisplay(
            recall=recall[i],
            precision=precision[i],
            average_precision=average_precision[i],
        )
        display.plot(ax=ax, name=f"Prec.-recall for {classes[i]}", color=color)
    ##### add the legend for the iso-f1 curves
    handles, labels = display.ax_.get_legend_handles_labels()
    handles.extend([l])
    labels.extend(["iso-f1 curves"])
    # set the legend and the axes
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.legend(handles=handles, labels=labels, loc="best")
    ax.set_title("Precision-Recall per relationship with \n score normalisation = {}".format(normalisation))
    plt.show()
    plt.close()
    plt.savefig(f'./data/precision_recall_per_relationship_class.png', format='png', dpi=300)

def plot_ROC(Y_test: np.array, y_score, normalisation = 'per sentence') -> None:
    """
    :param Y_test:
    :param y_score:
    :param normalisation: # legacy from previous code where normalisation method was varied
    :return: plots and saves ROC curves
    """
    n_classes = Y_test.shape[1]

    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(Y_test[:, i], y_score[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])

    # Compute micro-average ROC curve and ROC area
    fpr["micro"], tpr["micro"], _ = roc_curve(Y_test.ravel(), y_score.ravel())
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

    # First aggregate all false positive rates
    all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))

    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(n_classes):
        mean_tpr += np.interp(all_fpr, fpr[i], tpr[i])

    # Finally average it and compute AUC
    mean_tpr /= n_classes

    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

    # Plot all ROC curves
    plt.figure()
    plt.plot(
        fpr["micro"],
        tpr["micro"],
        label="micro-average ROC curve (area = {0:0.2f})".format(roc_auc["micro"]),
        color="deeppink",
        linestyle=":",
        linewidth=4,
    )

    plt.plot(
        fpr["macro"],
        tpr["macro"],
        label="macro-average ROC curve (area = {0:0.2f})".format(roc_auc["macro"]),
        color="navy",
        linestyle=":",
        linewidth=4,
    )

    colors = cycle(["navy", "turquoise", "darkorange", "cornflowerblue", "teal"])#cycle(["aqua", "darkorange", "cornflowerblue"])
    ##### I added classes name here
    classes = ['(N) inhibits', '(A+) agonism', '(B) binding', '(A-) antagonism']
    lw = 2
    for i, color in zip(range(n_classes), colors):
        plt.plot(
            fpr[i],
            tpr[i],
            color=color,
            lw=lw,
            label="ROC curve of {0} (area = {1:0.2f})".format(classes[i], roc_auc[i]),
        )

    plt.plot([0, 1], [0, 1], "k--", lw=lw)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC per relationship class with \n score normalisation = {}".format(normalisation))
    plt.legend(loc="lower right")
    plt.show()
    plt.close()
    plt.savefig(f'./data/ROC_per_relationship_class.png', format='png', dpi=300)

def main(nx_gnbr_merged: pd.DataFrame = None):
    """
    Runs all functions to plots curves
    :param nx_gnbr_merged:
    :returns: plots and saves figures
    """
    ##### load the gnbr_nexus merged data -this has 157646 rows and process for sci-kit learn #####
    if nx_gnbr_merged is not None:
        df = nx_gnbr_merged
    else:
        logger.info('loading merged data')
        df = pd.read_pickle(r'./data/df_gnbr_with_themes_not_normalised_merged_to_nexus_169k.pkl')  # this has 168076 rows

    df = explode_drop_duplicates(df)
    df = normalisation_fraction_across_all_themes_per_row(df)
    df_nx = mk_nx_ground_truth_labels(df)
    df_jn = join_themes_gnbr_to_nexus(df, df_nx)

    ##### convert df for sci-kit learn #####
    Y_test, y_score = get_yscore_ytest(df_jn)

    ##### Plotting section #####
    logger.info('plotting precision-recall and ROC')
    plot_counts_per_class_nx(df_jn)
    plot_precision_recall(Y_test, y_score)
    plot_ROC(Y_test, y_score)

if __name__ == "__main__":
    main()


