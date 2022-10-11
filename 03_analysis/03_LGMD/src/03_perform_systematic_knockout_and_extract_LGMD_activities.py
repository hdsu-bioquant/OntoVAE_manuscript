import os
import numpy as np
import pandas as pd
import pickle

from onto_vae.ontobj import *
from onto_vae.vae_model import *

project_dir = os.getcwd()

# import ontobj
with open(project_dir + '/ontobj/HPO_ontobj.pickle', 'rb') as f:
    hpo = pickle.load(f) 


# initialize OntoVAE model
hpo_model = OntoVAE(ontobj=hpo,
                dataset='recount3_GTEx',
                top_thresh=1000,
                bottom_thresh=10)
hpo_model.to(hpo_model.device)


# load best model
checkpoint = torch.load(project_dir + '/models/GTEx_HPO/best_model.pt', 
                        map_location=torch.device(hpo_model.device))
hpo_model.load_state_dict(checkpoint['model_state_dict'])


# extract genes
onto_genes = hpo.extract_genes(top_thresh=1000,
                               bottom_thresh=10)



# perform systematic knockout
ko_act = []

for g in onto_genes:

    act = hpo_model.perturbation(ontobj=hpo,
                                 dataset='recount3_GTEx',
                                 genes=[g],
                                 values=[0])
    ko_act.append(act)

ko_act = np.vstack(ko_act).T 


# load GTEx annotation
sample_annot = pd.read_csv(project_dir + '/datasets/recount3_GTEx/recount3_GTEx_annot.csv')
muscle_idx = sample_annot[sample_annot.tissue == 'Muscle'].index.to_numpy()


# find genes that were not 0 in all muscle samples in dataset
data = hpo.extract_data('recount3_GTEx',
                        top_thresh=1000,
                        bottom_thresh=10)
muscle_data = data[muscle_idx,:]
good_gene_ids = np.where(np.sum(muscle_data, axis=0) != 0)[0]

# subset
ko_act = pd.DataFrame(ko_act[muscle_idx, good_gene_ids])
ko_act.columns = np.array(onto_genes)[good_gene_ids]

# save results
ko_act.to_csv(project_dir + '/analysis/02_GTEx_HPO/results/systematic_knockout_GTEx_LGMD_activities.csv', index=False)


