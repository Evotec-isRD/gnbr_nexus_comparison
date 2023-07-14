#  A large-scale evaluation of NLP derived chemical-proteins relationships (GNBR vs Nexus): code for the paper

## Overview

This is code that accompanies a paper that is currently being submitted to PLOS ONE. 

The title of the paper is:
"A large-scale evaluation of NLP-derived chemical-gene/protein relationships from the scientific literature: 
implications for knowledge graph construction" 
(this will be a link to the paper DOI when it is published)

Our paper compares chemicals (small molecules) and proteins 'triples', such as "aspirin-INHIBITS-cyclooxygenase" between 
the NLP-derived dataset ["the Global Network of Biological Relationships (GNBR)"](https://academic.oup.com/bioinformatics/article/34/15/2614/4911883) 
and Evotec's own small molecule database called "Nexus" (how Nexus is constructed is outlined in our paper). 

## Description of code repository

The main.py file in the src directory executes all the other scripts. These scripts include preprocessing both GNBR 
and Nexus raw data for activities such as entity normalisation, then others for comparison and plotting 
(which script does what should be fairly obvious from their names). 

## Dependencies 
The library dependencies are in "requirements.txt"

By default, the code will run (executing from main.py) using data file dependencies that are a sub-set

There are several data files which are pre-requisites for the code to run. Where appropriate the external links 
to these files are provided. These need to be downloaded toa  "./data" directory that will have to be created locally. 
1. From the [GNBR database (version 7 2019)](https://zenodo.org/record/3459420#.Y8pCz3bP2Uk) the files 
**'part-ii-dependency-paths-chemical-gene-sorted-with-themes.txt.gz'** and **'part-i-chemical-gene-path-theme-distributions.txt.gz'**. 
These two files are mapped together in the gnbr_preprocess script. 
2. A **'gene_mapping_table.csv'** file, while is used to map Entrez gene IDs to UniProt IDs in GNBR data. 
3. A dictionary which has every chemical raw string name in gnbr 'fuzzy matched' to a Nexus chemical synonym. 
This file is called **'df_chemGeneSortedWithTheme_dropDup_NexusFuzzMatch_unzip.pkl'**. This is needed to 'normalise' gnbr names 
which are the same entity but have slightly different spellings from different authors in different papers - like an extra hyphen for instance. 
So an example in Nexus is 'protoporphyrin IX’ while in gnbr this could be 'protoporphyrin 9’. 
4. Nexus data in the file **'df_nexus_protein_level_averages.pkl'**. This has this has 2,462,326 rows and contains compound IDs 
linked with protein targets with an associated assay value (such as EC50). This is proprietary Evotec 
property and may be made available (in part) on request. 
5. A file that links Nexus chemical IDs to all known synonyms called **'structureID.csv'**
6. (this file is the result of 4 and 5 data merged together and processed with nexus_preprocess.py included here as a 
dependency to speed up processing time of the entire script). 
Nexus data processed to include a compound ID with all known name synonyms of that compound 
(if the compound has a name such as a drug), called **'pddf_nx_annotation_proteinLevelAverages_explode.parquet'**

## Instructions on running the code
1. Clone this project using your standard process
2. pip install the dependencies: "pip install -r requirements.txt" (use python=3.9)
3. run the main.py file: "python src/main.py" (this could take up to 30 minutes to run - watch the log!)

## Outputs
There should be several outputs on running the script:
1. Figures 1 then 3-5 and supporting information Figures. These will show via Matplotlib and also save to the 'data' folder in 
png format
2. Figure 2 was made manually from values outputted via the 'value_for_venn_figures.py' script, which appear in the console
via the logging library. 

## Contact
charlie.jeynes@evotec.com or tim.james@evotec.com

## Authors and acknowledgment
The code is by JCG (Charlie) Jeynes, and was reviewed by Matt Corney and Clare Eckold (both Evotec).
The paper is written by Charlie Jeynes, and was edited by Matt Corney and Tim James (both Evotec).

## License
The licence is a GNU General Public Licence (see the Licence.md file)

## Project status
This is in the process of being submitted to PLOS ONE. 
