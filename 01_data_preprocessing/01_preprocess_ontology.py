# This script shows how to preprocess a given ontology with the onto-vae package

# import all functions from onto-vae package
from onto_vae import *

# onto-vae takes as input two files:
# 1. a given ontology encoded as obo file
# 2. a tab-separated 2-column file containing the mapping from gene to GO term, the file would look like this
# ENSG00000177144	GO:0003723
# ENSG00000177144	GO:0005829
# ENSG00000177144	GO:0046872


# define a working directory
project_dir = '/user/ontologies/GO'
# onto-vae will create subfolders in the specified working directory and store its results there

# specify prefix to use in filenames for saving files
prefix='GO_ensembl'

# we initialize a new ontoobj 
# obo: path to our ontology obo file
# gene_annot: path to our 2-column gene -> term mapping file
# working_dir: where to store the results
# prefix: prefix used in filenames
onto = ontoobj(obo=onto_dir + '/data/go-basic.obo',
                      gene_annot=onto_dir + '/data/ensembl_goterm_mapping.txt',
                      working_dir=onto_dir,
                      prefix=prefix)

# The function 'create_dag_dict_files' will create the some initial files needed to train OntoVAE models
# id is an additional argument that can be passed if only a certain id from the ontology should be used
onto.create_dag_dict_files(id = 'biological_process')

# The function 'create_trim_dag_files' creates trimmed versions of the initial files using thresholds chosen by the user
onto.create_trim_dag_files(top_thresh=1000, bottom_thresh=30)

# The function 'create_model_input' creates a list of binary masks for the decoder part of OntoVAE
onto.create_model_input()

# The function 'compute_wsem_sim' creates a matrix of Wang semantic similarities (not necessarily needed by the user)
onto.compute_wsem_sim()