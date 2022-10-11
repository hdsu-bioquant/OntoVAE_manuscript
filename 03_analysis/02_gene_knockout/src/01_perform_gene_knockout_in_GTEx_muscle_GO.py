import os
import numpy as np
import pandas as pd
import pickle

from onto_vae.ontobj import *
from onto_vae.vae_model import *

project_dir = os.getcwd()


# import ontobj
with open(project_dir + '/ontobj/GO_ensembl_ontobj.pickle', 'rb') as f:
    go = pickle.load(f) 


# initialize OntoVAE model
go_model = OntoVAE(ontobj=go,
                dataset='recount3_GTEx')
go_model.to(go_model.device)


# load best model
checkpoint = torch.load(project_dir + '/models/GTEx_GO/best_model_neuronnum3_1.pt',
                        map_location=torch.device(go_model.device))
go_model.load_state_dict(checkpoint['model_state_dict'])


# extract genes 
onto_genes = go.extract_genes()


# load GTEx annotation
sample_annot = pd.read_csv(project_dir + '/datasets/recount3_GTEx/recount3_GTEx_annot.csv')
muscle_idx = sample_annot[sample_annot.tissue == 'Muscle'].index.to_numpy()

# extract dataset
data = go.extract_data('recount3_GTEx')
data = go.data['1000_30']['recount3_GTEx']

# subset to muscle samples
muscle_data = data[muscle_idx,:]

# add subsetted dataset to Ontoobj
go.add_dataset(dataset=data,
                description='recount3_GTEx_muscle',
                top_thresh=1000,
                bottom_thresh=10)
go.data['1000_30']['recount3_GTEx_muscle'] = muscle_data

# genes that the knockout will be performed on
ko_genes = ['ENSG00000198947', 'ENSG00000061936', 'ENSG00000178741'] # DMD, SFSWAP, COX5A

# helper function to perform knockout
def perform_ko(gene, output):
    ko_act = go_model.perturbation(ontobj=go,
                                 dataset='recount3_GTEx_muscle',
                                 genes=[gene],
                                 values=[0],
                                 output=output)
    return ko_act

# retrieve activities post knockout
DMD_act = perform_ko(ko_genes[0], 'terms')
np.save(project_dir + '/analysis/02_gene_knockout/results/activities/recount3_GTEx_muscle_DMD_knockout_GO_activities.npy', DMD_act)
SFSWAP_act = perform_ko(ko_genes[1], 'terms')
np.save(project_dir + '/analysis/02_gene_knockout/results/activities/recount3_GTEx_muscle_SFSWAP_knockout_GO_activities.npy', SFSWAP_act)
COX5A_act = perform_ko(ko_genes[2], 'terms')
np.save(project_dir + '/analysis/02_gene_knockout/results/activities/recount3_GTEx_muscle_COX5A_knockout_GO_activities.npy', COX5A_act)

# retrieve reconstructed values post knockout
DMD_rec = perform_ko(ko_genes[0], 'genes')
np.save(project_dir + '/analysis/02_gene_knockout/results/activities/recount3_GTEx_muscle_DMD_knockout_reconstructions.npy', DMD_rec)
SFSWAP_rec = perform_ko(ko_genes[1], 'genes')
np.save(project_dir + '/analysis/02_gene_knockout/results/activities/recount3_GTEx_muscle_SFSWAP_knockout_reconstructions.npy', SFSWAP_rec)
COX5A_rec = perform_ko(ko_genes[2], 'genes')
np.save(project_dir + '/analysis/02_gene_knockout/results/activities/recount3_GTEx_muscle_COX5A_knockout_reconstructions.npy', COX5A_rec)