import pandas as pd
import time
import itertools
from rapidfuzz import fuzz
from typing import Tuple
# import seaborn as sns
import matplotlib.pyplot as plt
import pandas_utilities as pu
from data_dependencies import DataDependencies
# sns.set(style='ticks', context='talk')
# sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.size': 22})
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

################################ unique raw string vs unique IDs (Figure 1) #####################################
def get_gnbr_unq_genesStr_vs_geneID(gnbr: pd.DataFrame) -> Tuple[int, int]:
    """
    unique gene names in GNBR, from the raw string name, and the entrezID
    how many gene name raw string have not been assigned a entrezID? 58300-47118 = 11182
    :returns gene_str (is 58300) , gene_ID (is 47118)
    """
    gene_str = len(gnbr['Second entity name, raw string'].astype(str).str.lower().unique())
    gene_ID = len(gnbr['Second entity name, database ID(s)'].astype(str).str.lower().unique())
    return len(gene_str), len(gene_ID)

def get_gnbr_unq_chemStr_vs_chemID(gnbr: pd.DataFrame) -> Tuple[int, int]:
    """
    unique gene names in GNBR, from the raw string name, and the entrezID
    how many gene name raw string have not been assigned a entrezID? 58300-47118 = 11182
    :returns gene_str (is 58300) , gene_ID (is 47118)
    """
    chem_str = gnbr['First entity name, raw string'].astype(str).str.lower().unique()
    chem_ID = gnbr['First entity name, database ID(s)'].astype(str).str.lower().unique()
    return len(chem_str), len(chem_ID)

def plot_raw_string_counts_to_IDs_for_paper(df: pd.DataFrame,
                                  cols1: list = ['First entity name, raw string','Second entity name, raw string'],
                                  cols2: list = ['First entity name, database ID(s)','Second entity name, database ID(s)'],
                                  title: str = None,
                                  ):
    """
    :param df:
    :param cols1: column names
    :param cols2: column names
    :param title:
    :return: plot and saved figure
    """
    df = pu.make_all_lower_case(df)
    df = pu.remove_leading_trailing_white_space(df)
    ##### Plot comparison of raw strings that are linked to an ID in GNBR
    group1data = df.loc[:, cols1].nunique().to_list()
    group1data = [64670, 58299]
    group1name = 'raw_string'
    group2data = df.loc[:, cols2].nunique().to_list()
    print(group2data)
    group2name = 'databaseID'
    ylabel = ['chemicals', 'genes/proteins']
    df_unqs = pd.DataFrame({group1name: group1data, group2name: group2data}, index=ylabel)
    ax = df_unqs[::-1].plot.barh()
    ax.set_xlabel('number of unique names')  # , fontdict={'fontsize': 30}
    ax.set_title(title)
    plt.tight_layout()
    plt.show()
    plt.close()
    plt.savefig(f'./data/raw_string_counts_to_IDs_for_paper.png', format='png', dpi=300)

####################################### chemicals (Figure 2) ################################################

class Chemicals:
    """get values for unique chemicals in gnbr and nexus and intersection between the two (Figure 2a)"""
    @staticmethod
    def get_nexus_number_unique_chemicals(nexus: pd.DataFrame) -> int:
        ''' with both the '/df_nexus_protein_level_averages.pkl' , and
        the '/df_nx_protein_level_averages_with_chemical_names.pkl',
        there are 2,462,326 'STRUCTURE_ID'. the unique of this is 1,151,388
        1113420 after lower case and removing trailing white space.
        Same number as the latter for and pddf_nx_annotation_proteinLevelAverages_explode.parquet'
        '''
        unq_chemicals = nexus['STRUCTURE_ID'].nunique()
        # formID = nexus['STRUCTURE_FORM_ID'].nunique() #same as above just checking
        logger.info(f'number of unique chemicals in nexus is {unq_chemicals}')
        return unq_chemicals

    @staticmethod
    def get_gnbr_number_unique_chemicals(gnbr: pd.DataFrame) -> int:
        """
        The number of unique chemical names in gnbr is 73273 (gnbr['First entity name, raw string'].unique()), but many of
        these are the same entity just with a hyphen or space for instance.
        This function matches names in gnbr to each other with a 'fuzzy-match' threshold of 0.95 (takes about 10 mins on
        a single CPU).
        This results in 8603 matches (i.e. these are the same entity)
        The out put of this function is 64,672, which is the value in Figure 2.
        """

        try:
            logger.info("loading ./data/s_gnbr_FuzzychemicalNamesSet.pkl")
            s_gnbr_FuzzychemicalNamesSet = pd.read_pickle(r'./data/s_gnbr_FuzzychemicalNamesSet.pkl')
        except:
            logger.info("making file ./data/s_gnbr_FuzzychemicalNamesSet.pkl. This could take up to 10 minutes")
            def find_similar_pairs(NamesLst, *, required_similarity=95):
                """
                Find pairs of similar-looking tags in the collection ``tags``.

                Increase ``required_similarity`` for stricter matching (=> less results).
                """
                for t1, t2 in itertools.combinations(sorted(NamesLst), 2):
                    if fuzz.ratio(t1, t2) > required_similarity:
                        yield [t1, t2]

            start = time.time()
            gnbr_FuzzychemicalNamesSet = list(
                find_similar_pairs(list(gnbr['First entity name, raw string'].astype(str).unique()),
                                   required_similarity=95)) # this takes about 10 mins
            print(time.time() - start)
            # convert to series
            s_gnbr_FuzzychemicalNamesSet = pd.Series(gnbr_FuzzychemicalNamesSet)
            # save it as a pickle file
            s_gnbr_FuzzychemicalNamesSet.to_pickle(r'./data/s_gnbr_FuzzychemicalNamesSet.pkl')
        logger.info(f"unique chemical names in gnbr is "
                    f"{len(gnbr['First entity name, raw string'].unique()) - len(s_gnbr_FuzzychemicalNamesSet)} ")

        return len(gnbr['First entity name, raw string'].unique()) - len(s_gnbr_FuzzychemicalNamesSet)

    @staticmethod
    def intersection_chem_names_nx_gnbr()-> int:
        """
        The intersection of chemical names in nexus and GNBR. This is calculated by counting the length of
        fuzzy matched pairs between the two databases.
        # length just synonyms = 18373
        # length including mechaniams of action = 22426
        # length including mechaniams of action + molecule type = 22426
        """
        # get fuzzy match dictionary between nexus and gnbr chem names
        fuzzMatch_nexus_gnbrWithThemesDrpDup = pd.read_pickle(r'./data/df_chemGeneSortedWithTheme_dropDup_NexusFuzzMatch_unzip.pkl')
        # extract 1st tuple (chemical name) from what rapid fuzz returns
        fuzzMatch_nexus_gnbrWithThemesDrpDup['1st_tuple'] = fuzzMatch_nexus_gnbrWithThemesDrpDup.iloc[:, 2].apply(
            lambda x: x if x is None else x[0])
        # counts the number of unique chemical names that have been matched to gnbr
        chem_intersection = fuzzMatch_nexus_gnbrWithThemesDrpDup[
            fuzzMatch_nexus_gnbrWithThemesDrpDup['1st_tuple'].notnull()]
        chem_intersection = chem_intersection['1st_tuple'].astype(str).nunique()
        logger.info(f'the intersection in chemical names between gnbr and nexus is {chem_intersection}')
        return chem_intersection

#################################### gene/proteins (Figure 2) #############################################
class GenesProteins:
    """get values for unique chemicals in gnbr and nexus and intersection between the two (Figure 2a)"""

    @staticmethod
    def get_unique_proteins_nexus(nexus: pd.DataFrame):
        """Get the unique proteins in nexus using the UniProt ID. This is 14,698 unique IDs. This function also
        gets the percentage of proteins by species which is
        Homo sapiens                      0.817241
        Rattus norvegicus                 0.061071
        Mus musculus                      0.022868
        Photinus pyralis                  0.008925
        """
        # get unique protein IDs
        unq_proteins = len(nexus['UNIPROT_ACCESSION_NUMBERS'].unique())
        # get the percentage of proteins that are from different species
        percentage_species = nexus['SPECIES'].value_counts(normalize=True)

        return unq_proteins

    @staticmethod
    def get_unique_proteins_gnbr(gnbr: pd.DataFrame):
        """Get the unique proteins in gnbr using the UniProt ID. This is 47118 unique IDs. This function also
        gets the percentage of proteins by species which is
        human         0.642549
        tax:10116)    0.156372
        tax:10090)    0.137725
        tax:9913)     0.008912
        tax:4932)     0.008435
        """
        # get unique protein IDs
        unq_proteins = len(gnbr['Second entity name, database ID(s)'].astype(str).str.lower().unique())
        # get the percentage of proteins that are from different species
        # percentage_species  = gnbr['taxonomy'].value_counts(normalize=True) # works but need the columns not filtered
        logger.info(f'the number of unique genes/ protein in gnbr is {unq_proteins}')
        return unq_proteins

    @staticmethod
    def get_intersecting_proteins_gnbr_nexus(gnbr: pd.DataFrame, nexus: pd.DataFrame) -> int:
        """
        To do this first the gnbr gene IDs from Entrez have to be converted to UniProt IDs. This is done
        via a mapping file that is loaded in. Then the UniProt IDs can be compared to those in nexus.
        :return: the intersection is 3764 proteins
        """
        # Load geneMapping file
        df_geneMapping = pd.read_csv(r'./data/gene_mapping_table.csv')
        # convert P_ENTREZGENEID ints to str
        df_geneMapping['P_ENTREZGENEID'] = df_geneMapping['P_ENTREZGENEID'].astype('Int64').astype(str)
        # make sure the uniprot ID are lower cases as all gnbr is lower case
        # map
        dict_entrezID_to_uniprot = dict(zip(df_geneMapping['P_ENTREZGENEID'], df_geneMapping['ACC']))
        # create uniprot column in gnbr
        gnbr['uniprotID'] = gnbr['Second entity name, database ID(s)'].map(dict_entrezID_to_uniprot)
        # make uniprot column lower case in gnbr because we did lower case earlier to nexus
        gnbr['uniprotID'] = gnbr['uniprotID'].str.lower()
        # find intersection
        nexus_uniprot_set = set([x for x in nexus['UNIPROT_ACCESSION_NUMBERS']])
        gnbr_uniprot_set = set([x for x in gnbr['uniprotID']])
        return len(nexus_uniprot_set.intersection(gnbr_uniprot_set))

############################## chemical-gene/proteins pairs (Figure 2c) #######################################

class ChemicalGenePairs:
    """get values for chemical-gene/proteins pairs (Figure 2c)"""
    @staticmethod
    def get_unique_ChemProteinsPairs_gnbr(gnbr: pd.DataFrame) -> int:
        """Get the unique protein-chemical pairs from gnbr.
        unique chem_MSH_ID to gene_ID pairs: this is 158,625
        """
        # note: the following lines experimented with normalising chemical names with fuzz matching
        # get fuzz match names that normalise on the same entities which may have things like an extra hyphen
        # fuzzymatched_gnbr_chem_names = pd.read_csv(r'./data/s_gnbr_FuzzychemicalNamesSet.pkl')
        # need to have a good way to map correctly
        # map fuzzy matched names to normalise on the same entities which may have things like an extra hyphen
        # gnbr['fuzzed_chem_name'] = gnbr['First entity name, raw string'].map()

        # get it chem_RAW_STR to gene_ID: this is 370,897
        unique_ChemstrProteinsPairs_gnbr = (gnbr.
                                                drop_duplicates(
            subset=['First entity name, raw string',  'Second entity name, database ID(s)'], keep='last')
                                            .dropna()
                                            )


        # get it chem_MSH_ID to gene_ID: this is 326,066
        unique_ChemIDProteinsPairs_gnbr = (gnbr.
                                                drop_duplicates(
            subset=['First entity name, database ID(s)',  'Second entity name, database ID(s)'], keep='last')
                                            .dropna()
                                            )
        logger.info(f'number of unique chem-protein pairs in gnbr is {len(unique_ChemIDProteinsPairs_gnbr)}')

        return len(unique_ChemIDProteinsPairs_gnbr)

    @staticmethod
    def get_unique_ChemProteinsPairs_nexus(nexus: pd.DataFrame) -> int:
        """Get the unique protein-chemical pairs from nexus.
        unique structure_ID to uniprot_ID pairs: this is 2,738,658 (2,766,995 not str lowered)
        """

        # unique structure_ID to uniprot_ID pairs: this is 2,766,995
        unique_ChemIDProteinsPairs_nexus = (nexus.
                                                drop_duplicates(
            subset=['STRUCTURE_ID', 'UNIPROT_ACCESSION_NUMBERS',], keep='last')
                                             )

        logger.info(f'number of unique chem-protein pairs in nexus is {len(unique_ChemIDProteinsPairs_nexus)}')


        return len(unique_ChemIDProteinsPairs_nexus)

    @staticmethod
    def get_unique_intersecting_ChemProteinsPairs() -> int:
        """Get the unique protein-chemical pairs that are present in both nexus and gnbr.
        This loads in a file, from a merge done on fuzzy name matches to both gnbr and nexus
        and then gets the unique chemical names paired with a UniProt ID. Note this file is created in
        the merge_nexus_gnbr.py file.
        the nexus cols give 9,990 and the gnbr cols give 10,699
        take the nexus as they have unq structure ID whereas the gnbr is likely to have 'two' entities that are actually
        the same
        """
        df_mg = pd.read_pickle('./data/df_gnbr_with_themes_not_normalised_merged_to_nexus_169k.pkl')
        # nexus cols gives 9,990 (9991 with 'nan' included)
        unq_intersecting_chem_gene_pairs_nx_cols = df_mg.drop_duplicates(subset=['chem_ID', 'protein_ID'])

        logger.info(f'number of intersecting unique chem-protein pairs in gnbr and nexus is {len(unq_intersecting_chem_gene_pairs_nx_cols)}')

        # gnbr cols 10,699 this is just for comparison
        unq_intersecting_chem_gene_pairs_gnbr_cols = len(
            df_mg.drop_duplicates(subset=['chemical, raw string', 'uniprotID']))


        return len(unq_intersecting_chem_gene_pairs_nx_cols)

#################################### loading data #############################################
def get_nexus_data_with_exploded_annotations():
    """
     Load in nexus data with proteins levels joined to annoatations
     has 19,406,072 exploded rows (11,536,091 rows before explosion)
    """
    pddf_nx_merge_explode = pd.read_parquet(DataDependencies.NEXUS_PROCESSED)

    pddf_nx_merge_explode = pddf_nx_merge_explode[['STRUCTURE_ID', 'UNIPROT_ACCESSION_NUMBERS', 'TARGET_NAME', \
                                                   'TARGET_TYPE', 'SPECIES', 'MOLECULE_FUNCTION', 'ID',
                                                   'STRUCTURE_FORM_ID', \
                                                   'DESCRIPTION', 'ANNOTATION', 'CREATED']]

    return pddf_nx_merge_explode

def get_gnbr():
    """
    Loads data and selects columns. Gives 1728361 rows.
    """
    gnbr = pd.read_csv(r'./data/part-ii-dependency-paths-chemical-gene-sorted-with-themes.txt.gz', header=None, sep='\t')

    column_names = ['PubMed ID',
                    'Sentence number (0 = title)',
                    'First entity name, formatted',
                    'First entity name, location (characters from start of abstract)',
                    'Second entity name, formatted',
                    'Second entity name, location',
                    'First entity name, raw string',
                    # 'chemical, raw string',
                    'Second entity name, raw string',
                    # 'gene, raw string',
                    'First entity name, database ID(s)',
                    # 'chemical, database ID(s)',
                    'Second entity name, database ID(s)',
                    # 'gene, database ID(s)',
                    'First entity type (Chemical, Gene, Disease)',
                    'Second entity type (Chemical, Gene, Disease)',
                    'Dependency path',
                    'Sentence, tokenized']

    gnbr.columns = column_names

    # select only a few columns
    gnbr = gnbr[
        ['PubMed ID',
         'First entity name, raw string',
         'Second entity name, raw string',
         'First entity name, database ID(s)',
         'Second entity name, database ID(s)',
         'First entity type (Chemical, Gene, Disease)',
         'Second entity type (Chemical, Gene, Disease)',
         'Dependency path',
         'Sentence, tokenized']]


    return gnbr

def main():
    # load data
    gnbr = get_gnbr()
    nexus = get_nexus_data_with_exploded_annotations()

    # drop cases etc
    gnbr = pu.make_all_lower_case(gnbr)
    gnbr = pu.remove_leading_trailing_white_space(gnbr)
    nexus = pu.make_all_lower_case(nexus)
    nexus = pu.remove_leading_trailing_white_space(nexus)

    # Figure 1
    plot_raw_string_counts_to_IDs_for_paper(gnbr)

    # values for Figure 2 Venn plot
    Chemicals.get_gnbr_number_unique_chemicals(gnbr) # 64,672 gnbr chems
    Chemicals.get_nexus_number_unique_chemicals(nexus) # 1,113,420 nexus chems
    Chemicals.intersection_chem_names_nx_gnbr() # intersecting chems 18,045

    GenesProteins.get_unique_proteins_gnbr(gnbr) # unique gnbr genes/proteins = 47,118
    GenesProteins.get_unique_proteins_nexus(nexus) # unique nexus genes/proteins = 14,698
    GenesProteins.get_intersecting_proteins_gnbr_nexus(gnbr, nexus) # intersecting gene/proteins = 4,124

    ChemicalGenePairs.get_unique_ChemProteinsPairs_gnbr(gnbr) # unique gnbr genes/proteins-chem pairs = 158,625
    ChemicalGenePairs.get_unique_ChemProteinsPairs_nexus(nexus) # unique nexus genes/proteins-chem pairs = 2,738,658
    ChemicalGenePairs.get_unique_intersecting_ChemProteinsPairs() # unique nexus genes/proteins-chem pairs = 9,990

#
if __name__ == '__main__':
    main()
