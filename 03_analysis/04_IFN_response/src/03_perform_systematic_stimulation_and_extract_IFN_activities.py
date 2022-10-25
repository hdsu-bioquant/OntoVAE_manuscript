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



# load sample annotation
sample_annot =  pd.read_csv(project_dir + '/datasets/Kang_PBMC/Kang_PBMC_control_annot.csv')
cd4t_idx = sample_annot[sample_annot.celltype == 'CD4T'].index.to_numpy()

# extract dataset
data = go.extract_dataset('Kang_PBMC_control',
                            top_thresh=1000,
                            bottom_thresh=30)


# subset to CD4T cells
cd4t_data = data[cd4t_idx,:]

# get genes whose expression is not 0 for all cells
good_gene_ids = np.where(np.sum(data, axis=0) != 0)[0]
good_genes = np.array(onto_genes)[good_gene_ids]

# add subsetted dataset to Ontoobj
go.add_dataset(dataset = cd4t_data,
               description = 'CD4T_PBMC_control',
               top_thresh=1000,
               bottom_thresh=30)


# perform systematic stimulation of good genes in CD4T cells and 
stim_act = np.zeros((cd4t_data.shape[0],len(good_genes)))

for i in range(len(good_genes)):   
    print(good_genes[i])
    act = go_model.perturbation(ontobj=go,
                                 dataset='CD4T_PBMC_control',
                                 genes=[good_genes[i]],
                                 values=[8],
                                 terms=['GO:0060337'])
    stim_act[:,i] = act.copy().flatten()


# convert to dataframe
stim_act = pd.DataFrame(stim_act)
stim_act.columns = np.array(good_genes)

# save results
stim_act.to_csv(project_dir + '/03_analysis/04_IFN_response/results/systematic_stimulation_control_CD4T_IFN_activities.csv', index=False)