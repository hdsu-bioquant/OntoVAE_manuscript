# for HPO + recount3 GTEx

import os
import pickle
from onto_vae.ontobj import *

project_dir = os.getcwd()

# initialize the Ontobj
hpo = Ontobj(description='HPO')

# initialize ontology
hpo.initialize_dag(obo=project_dir + '/ontologies/HPO/hp.obo',
                   gene_annot=project_dir + '/ontologies/HPO/gene_term_mapping.txt')

# trim the ontology
hpo.trim_dag(top_thresh=1000, 
             bottom_thresh=10)

# create masks for decoder initialization
hpo.create_masks(top_thresh=1000,
                 bottom_thresh=10)

# match recount3 GTEx dataset to the ontology
hpo.match_dataset(expr_data = project_dir + '/datasets/recount3_GTEx/recount3_GTEx_symbol_log_tpm.csv',
                  name='recount3_GTEx',
                  top_thresh=1000,
                  bottom_thresh=10)

# save the ontobj
with open(project_dir + '/ontobj/HPO_ontobj.pickle', 'wb') as f:
    pickle.dump(hpo, f) 
