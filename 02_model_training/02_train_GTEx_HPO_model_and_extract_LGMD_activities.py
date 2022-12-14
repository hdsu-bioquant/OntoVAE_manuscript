import os
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

# train model
hpo_model.train_model(project_dir + '/models/GTEx_HPO/best_model.pt',
                     log=False)


# load best model
checkpoint = torch.load(project_dir + '/models/GTEx_HPO/best_model.pt', 
                        map_location=torch.device(hpo_model.device))
hpo_model.load_state_dict(checkpoint['model_state_dict'])


# retrieve pathway activities for LGMD
act = hpo_model.get_pathway_activities(ontobj=hpo,
                                       dataset='recount3_GTEx',
                                       terms=['HP:0006785']) 

# load GTEx annotation
sample_annot = pd.read_csv(project_dir + '/datasets/recount3_GTEx/recount3_GTEx_annot.csv')
muscle_idx = sample_annot[sample_annot.tissue == 'Muscle'].index.to_numpy()

# subset to muscle samples
act = act[muscle_idx,:]

# save LGMD activities
np.save(project_dir + '/activities/GTEx_muscle_LGMD_activities.npy', act.flatten())

