###-------------------------------------------------------------###
##        TRAINING OF ONTOLOGY VARIATIONAL AUTOENCODER           ##
###-------------------------------------------------------------###


import sys
import os
import re
import json
import pandas as pd 
import numpy as np
import torch 
from torch import optim
import torch.nn as nn
import torch.nn.functional as F
import neptune.new as neptune
import pickle
from tqdm import tqdm

from onto_vae import *


DATA = sys.argv[1] 
PARAMETERS = sys.argv[2] 
MASKS = sys.argv[3]
MODELPATH = sys.argv[4]
NAME = sys.argv[5]



###------------------------------------------
## DATA IMPORT AND PREPARATION
###------------------------------------------

# load counts
X = np.load(DATA)
samp_num, feat_num = X.shape

###------------------------------------------
## PARAMETER DEFINITION
###------------------------------------------

with open(PARAMETERS) as file:
    PARAMS = json.load(file)



###------------------------------------------
## DECODER INITIALIZATION
###------------------------------------------

# load binary matrix list for decoder connections
with open(MASKS, 'rb') as f:
    matrix_list = pickle.load(f)

# convert binary matrix list to tensor
mask_list = [torch.tensor(m, dtype=torch.float32) for m in matrix_list]

# define layer dims of decoder
layer_dims_onto = np.array([mask_list[0].shape[1]] + [m.shape[0] for m in mask_list])

###------------------------------------------
## MODEL CREATION
###------------------------------------------

# create model
model = OntoVAE(onto = PARAMS['onto'],
               in_features = feat_num,
               layer_dims_dec = layer_dims_onto,
               mask_list = mask_list,
               neuronnum = PARAMS['neuronnum'],
               drop = PARAMS['dropout'],
               z_drop = PARAMS['z_dropout'])

# move model to cuda
device = torch.device('cuda')
model.to(device)



###------------------------------------------
## LOAD BEST MODEL
###------------------------------------------

checkpoint = torch.load(MODELPATH + '/best_model.pt')
model.load_state_dict(checkpoint['model_state_dict'])



###------------------------------------------
## LATENT SPACE EMBEDDING
###------------------------------------------

X = torch.tensor(X, dtype=torch.float32).to(device)

model.eval()
with torch.no_grad():
    z = model.get_embedding(X)
    z = z.to('cpu').detach().numpy()

z = np.array(np.split(z, z.shape[1]/PARAMS['neuronnum'], axis=1)).mean(axis=2).T
np.save(MODELPATH + '/' + NAME + '_latent_space_embedding.npy', z)




###------------------------------------------
## PATHWAY ACTIVATIONS
###------------------------------------------

# hook function to extract activations

activation = {}
def get_activation(index):
    def hook(model, input, output):
        activation[index] = output.to('cpu').detach()
    return hook

# register forward hooks
hooks = {}

for i in range(len(model.decoder.decoder)-1):
    key = 'Dec' + str(i)
    value = model.decoder.decoder[i][0].register_forward_hook(get_activation(i))
    hooks[key] = value

# run data through model
model.eval()
with torch.no_grad():
    reconstruction, _, _ = model(X)

# concatenate activations
act = torch.cat(list(activation.values()), dim=1).detach().numpy()
act = np.array(np.split(act, act.shape[1]/PARAMS['neuronnum'], axis=1)).mean(axis=2).T

# save results
np.save(MODELPATH + '/' + NAME + '_decoder_activities.npy', act)
np.save(MODELPATH + '/' + NAME + '_reconstruction.npy', reconstruction.to('cpu').detach().numpy())
