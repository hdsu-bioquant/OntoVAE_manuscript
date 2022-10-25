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


# load GTEx annotation
sample_annot = pd.read_csv(project_dir + '/datasets/recount3_GTEx/recount3_GTEx_annot.csv')
muscle_idx = sample_annot[sample_annot.tissue == 'Muscle'].index.to_numpy()

# extract dataset
data = hpo.extract_data('recount3_GTEx',
                        top_thresh=1000,
                        bottom_thresh=10)

# subset to muscle samples
data = data[muscle_idx,:]

# get genes whose expression is not 0 for all samples
good_gene_ids = np.where(np.sum(muscle_data, axis=0) != 0)[0]
good_genes = np.array(onto_genes)[good_gene_ids]

# add subsetted dataset to Ontoobj
hpo.add_dataset(dataset=data,
                description='recount3_GTEx_muscle',
                top_thresh=1000,
                bottom_thresh=10)

# perform systematic knockout in muscle
ko_act = np.zeros((data.shape[0],len(good_genes)))

for i in range(len(good_genes)):

    act = hpo_model.perturbation(ontobj=hpo,
                                 dataset='recount3_GTEx_muscle',
                                 genes=[good_genes[i]],
                                 values=[0],
                                 terms=['HP:0006785'])
    ko_act[:,i] = act.flatten()


# create dataframe
ko_act = pd.DataFrame(ko_act)
ko_act.columns = good_genes

# save results
ko_act.to_csv(project_dir + '/analysis/03_LGMD/results/systematic_knockout_GTEx_LGMD_activities.csv', index=False)


