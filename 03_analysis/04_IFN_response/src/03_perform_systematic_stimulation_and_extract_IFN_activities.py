import os
import numpy as np
import pandas as pd
import pickle

from onto_vae.ontobj import *
from onto_vae.vae_model import *

project_dir = os.getcwd()

# import ontobj
with open(project_dir + '/ontobj/GO_symbol_ontobj.pickle', 'rb') as f:
    go = pickle.load(f) 


# initialize OntoVAE model
go_model = OntoVAE(ontobj=go,
                dataset='Kang_PBMC_control',
                top_thresh=1000,
                bottom_thresh=30)
go_model.to(go_model.device)


# load best model
checkpoint = torch.load(project_dir + '/models/PBMC_GO/best_model.pt', 
                        map_location=torch.device(go_model.device))
go_model.load_state_dict(checkpoint['model_state_dict'])


# extract genes
onto_genes = go.extract_genes(top_thresh=1000,
                               bottom_thresh=30)



# perform systematic stimulation
stim_act = []

for g in onto_genes:

    act = go_model.perturbation(ontobj=go,
                                 dataset='Kang_PBMC_control',
                                 genes=[g],
                                 values=[8])
    stim_act.append(act)

stim_act = np.vstack(stim_act).T 


# load Kang PBMC annotation
sample_annot = pd.read_csv(project_dir + '/datasets/Kang_PBMC/Kang_PBMC_control_annot.csv')
cd4t_idx = sample_annot[sample_annot.celltype == 'CD4T'].index.to_numpy()


# find genes that were not 0 in all muscle samples in dataset
data = go.extract_dataset('Kang_PBMC_control',
                        top_thresh=1000,
                        bottom_thresh=30)
cd4t_data = data[cd4t_idx,:]
good_gene_ids = np.where(np.sum(cd4t_data, axis=0) != 0)[0]

# subset
stim_act = pd.DataFrame(stim_act[cd4t_idx, good_gene_ids])
stim_act.columns = np.array(onto_genes)[good_gene_ids]

# save results
stim_act.to_csv(project_dir + '/analysis/04_IFN_response/results/systematic_stimulation_PBMC_CD4T_IFN_activities.csv', index=False)


act_frame_sub = act_frame.iloc[cd4t_idx,:]