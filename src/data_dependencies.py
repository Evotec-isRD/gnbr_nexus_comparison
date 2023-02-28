"""
There are several data files which are pre-requisites for the code to run. Please contact Evotec (see contact section in readme)
to get these files. These need to be in a "./data" directory that will have to be created locally.
1a. From the [GNBR database (version 7 2019)](https://zenodo.org/record/3459420#.Y8pCz3bP2Uk) the files
**'part-ii-dependency-paths-chemical-gene-sorted-with-themes.txt.gz'**
1b. and **'part-i-chemical-gene-path-theme-distributions.txt.gz'**. These two files are mapped together in the gnbr_preprocess script.
2. A **'gene_mapping_table.csv'** file, while is used to map Entrez gene IDs to UniProt IDs in GNBR data.
3. A dictionary which has every chemical raw string name in gnbr, 'fuzzy matched' to a Nexus chemical synonym.
This file is called **'df_chemGeneSortedWithTheme_dropDup_NexusFuzzMatch_unzip.pkl'**. This is needed to 'normalise' gnbr names
which are the same entity but have slightly different spellings from different authors in different papers - like an extra hyphen for instance.
So an example in Nexus is 'protoporphyrin IX’ while in gnbr this could be 'protoporphyrin 9’.
4. Nexus data in the file **'df_nexus_protein_level_averages.pkl'**. This has this has 2,462,326 rows and contains compound IDs
linked with protein targets with an associated assay value (such as EC50). This is proprietary Evotec
property and may be made available (in part) on request.
5. A file that links Nexus chemical IDs to all known synonyms called **'structureID.csv'**
6. (this file is the result of 4 and 5 data merged together and processed with nexus_preprocess.py
- included here as a dependency to speed up processing time of the entire script)
Nexus data processed to include a compound ID with all known name synonyms of that compound
(if the compound has a name such as a drug), called **'pddf_nx_annotation_proteinLevelAverages_explode.parquet'**
"""

class DataDependencies:
    """this contains all data files (paths to 'data' folder) needed to run the main.py script"""
    GNBR_DISTRIBUTIONS = r'./data/part-i-chemical-gene-path-theme-distributions.txt.gz' #  see 1a. in this script docstring
    GNBR_CHEMGENES = r'./data/part-ii-dependency-paths-chemical-gene-sorted-with-themes.txt.gz' #  see 1b. in this script docstring
    GENE_MAPPING = r'./data/gene_mapping_table.csv' #  see 2. in this script docstring
    FUZZMATCHED_CHEMNAME_NEXUS_GNBR = r'./data/df_chemGeneSortedWithTheme_dropDup_NexusFuzzMatch_unzip.pkl'  # see 3. in this script docstring
    NEXUS_UNPROCESSED = r'./data/df_nexus_protein_level_averages.pkl' # see 4. in this script docstring
    NEXUS_SYNONYMS = r'./data/structureID.csv' # see 5. in this script docstring
    NEXUS_PROCESSED = r'./data/pddf_nx_annotation_proteinLevelAverages_explode.parquet' # see 6. in this script docstring


