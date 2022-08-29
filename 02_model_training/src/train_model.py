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
MODELDIR = sys.argv[4]
LOG = sys.argv[5].lower() == 'true'


###------------------------------------------
## DATA IMPORT AND PREPARATION
###------------------------------------------

# load counts
X = np.load(DATA)
samp_num, feat_num = X.shape

# split data into train and val
indices = np.random.RandomState(seed=42).permutation(X.shape[0])
X_train_ind = indices[:round(len(indices)*0.8)]
X_val_ind = indices[round(len(indices)*0.8):]
X_train, X_val = X[X_train_ind,:], X[X_val_ind,:]

# convert train and val into torch tensors
X_train = torch.tensor(X_train, dtype=torch.float32)
X_val = torch.tensor(X_val, dtype=torch.float32)



###------------------------------------------
## PARAMETER DEFINITION
###------------------------------------------

with open(PARAMETERS) as file:
    PARAMS = json.load(file)



###------------------------------------------
## DATA LOADER
###------------------------------------------

X_train_batches = FastTensorDataLoader(X_train, 
                                       batch_size=PARAMS['batch_size'], 
                                       shuffle=True)
X_val_batches = FastTensorDataLoader(X_val, 
                                     batch_size=PARAMS['batch_size'], 
                                     shuffle=False)



###------------------------------------------
## LOGGER
###------------------------------------------

if LOG:
    run = neptune.init(
        project="neptune-project-here",
        api_token="api-token-here",
    ) 

    run["parameters"] = PARAMS
    run['sample_num'] = samp_num
    run['feat_num'] = feat_num

    MODELPATH = MODELDIR + str(vars(run)['_short_id']) 

    if not os.path.exists(MODELPATH):
        os.mkdir(MODELPATH) 

else:
    MODELPATH = MODELDIR

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
## MODEL TRAINING
###------------------------------------------

model.train_model(
    trainloader = X_train_batches,
    valloader = X_val_batches,
    lr = PARAMS['learning_rate'],
    kl_coeff = PARAMS['kl_coeff'],
    epochs = PARAMS['max_epochs'],
    modelpath = MODELPATH + '/best_model.pt',
    log = LOG,
    run = run
)


run.stop()
