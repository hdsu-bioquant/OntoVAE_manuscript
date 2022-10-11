# for GO + recount3 GTEx

import os
import pickle
from onto_vae.ontobj import *

project_dir = os.getcwd()

# initialize the Ontobj
go = Ontobj(description='GO_ensembl')

# initialize ontology
go.initialize_dag(obo=project_dir + '/ontologies/GO/go-basic.obo',
                   gene_annot=project_dir + '/ontologies/GO/ensembl_goterm_mapping.txt',
                   id = 'biological_process')

# trim the ontology
go.trim_dag(top_thresh=1000, 
             bottom_thresh=30)

# create masks for decoder initialization
go.create_masks(top_thresh=1000,
                 bottom_thresh=30)

# compute Wang Semantic similarities (used for app)
go.compute_wsem_sim(obo=project_dir + '/ontologies/GO/go-basic.obo',
                     top_thresh=1000,
                     bottom_thresh=30)

# match recount3 GTEx dataset to the ontology
go.match_dataset(expr_data = project_dir + '/datasets/recount3_GTEx/recount3_GTEx_ensembl_log_tpm.csv',
                  name='recount3_GTEx')

# save the ontobj
with open(project_dir + '/ontobj/GO_ensembl_ontobj.pickle', 'wb') as f:
    pickle.dump(go, f) 