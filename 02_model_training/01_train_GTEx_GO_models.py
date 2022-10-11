import os
import pickle
from onto_vae.ontobj import *
from onto_vae.vae_model import *

project_dir = os.getcwd()



# import ontobj
with open(project_dir + '/ontobj/GO_ensembl_ontobj.pickle', 'rb') as f:
    go = pickle.load(f) 

# helper function to train the model
def train_GTEx_GO_model(neuronnum, repeat):


    # initialize OntoVAE model
    go_model = OntoVAE(ontobj=go,
                    dataset='recount3_GTEx',
                    top_thresh=1000,
                    bottom_thresh=30,
                    neuronnum=neuronnum)
    go_model.to(go_model.device)

    # train model
    go_model.train_model(project_dir + '/models/GTEx_GO/best_model_neuronnum' + str(neuronnum) + '_' + repeat + '.pt',
                        log=False)


    # load best model
    checkpoint = torch.load(project_dir + '/models/GTEx_GO/best_model_neuronnum' + str(neuronnum) + '_' + repeat + '.pt', 
                            map_location=torch.device(go_model.device))
    go_model.load_state_dict(checkpoint['model_state_dict'])


    # retrieve pathway activities 
    act = go_model.get_pathway_activities(ontobj=go,
                                        dataset='recount3_GTEx')


    # save pathway activities
    np.save(project_dir + '/activities/recount3_GTEx_GO_activities_neuronnum' + str(neuronnum) + '_' + repeat + '.npy', act)



# train the model for neuronnums 1,2 and 3 twice
train_GTEx_GO_model(1, '1')
train_GTEx_GO_model(1, '2')
train_GTEx_GO_model(2, '1')
train_GTEx_GO_model(2, '2')
train_GTEx_GO_model(3, '1')
train_GTEx_GO_model(3, '2')