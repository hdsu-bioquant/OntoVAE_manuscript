import os
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

# train model
go_model.train_model(project_dir + '/models/PBMC_GO/best_model.pt',
                     log=False)


# load best model
checkpoint = torch.load(project_dir + '/models/PBMC_GO/best_model.pt',
                        map_location=torch.device(go_model.device))
go_model.load_state_dict(checkpoint['model_state_dict'])


# retrieve pathway activities for IFN signaling
act = go_model.get_pathway_activities(ontobj=go,
                                       dataset='Kang_PBMC_control',
                                       terms=['GO:0060337']) 

# load PBMC control annotation
sample_annot = pd.read_csv(project_dir + '/datasets/Kang_PBMC/Kang_PBMC_control_annot.csv')
cd4t_idx = sample_annot[sample_annot.celltype == 'CD4T'].index.to_numpy()

# subset to CD4T cells
act = act[cd4t_idx,:]

# save LGMD activities
np.save(project_dir + '/activities/PBMC_control_CD4T_IFN_activities.npy', act.flatten())
