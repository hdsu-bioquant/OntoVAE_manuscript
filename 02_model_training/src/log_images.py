###-------------------------------------------------------------###
##                    LOGGING IMAGES TO NEPTUNE                  ##
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
import umap
from sklearn.decomposition import PCA
import neptune.new as neptune
import pickle
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns 
import colorcet as cc

from onto_vae import *

DATA = sys.argv[1] 
PARAMETERS = sys.argv[2] 
MASKS = sys.argv[3]
RUN = sys.argv[4]
MODELPATH = sys.argv[5]
ONTO_ANNOT = sys.argv[6]
SAMPLE_ANNOT = sys.argv[7]
ANNOT_COLS = str.split(sys.argv[8], ',')
NAME = sys.argv[9]


###------------------------------------------
## LOGGER
###------------------------------------------

run = neptune.init(
    project="",
    run=RUN,
    api_token=""
) 




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

###------------------------------------------
## UMAP
###------------------------------------------

reducer = umap.UMAP()
embedding = reducer.fit_transform(z)

z = np.array(np.split(z, z.shape[1]/PARAMS['neuronnum'], axis=1)).mean(axis=2).T
np.save(MODELPATH + '/' + NAME + '_latent_space_embedding.npy', z)



###------------------------------------------
## UMAP PLOTTING
###------------------------------------------

fig, ax = plt.subplots(1,2, figsize=(20,10))
sns.scatterplot(x=embedding[:,0],
            y=embedding[:,1], 
            hue=annot.loc[:,ANNOT_COLS[0]],
            palette=color_dict,
            legend='full',
            s=8,
            ax=ax[0])
ax[0].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
ax[0].set_xlabel('UMAP1')
ax[0].set_ylabel('UMAP2')

if(len(ANNOT_COLS) > 1):
    sns.scatterplot(x=embedding[:,0],
                y=embedding[:,1], 
                hue=annot.loc[:,ANNOT_COLS[1]],
                palette=color_dict2,
                legend='full',
                s=8,
                ax=ax[1])
    ax[1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax[1].set_xlabel('UMAP1')
    ax[1].set_ylabel('UMAP2')

plt.tight_layout()
run["images/UMAP"].upload(fig)



###------------------------------------------
## LATENT SPACE DIM PLOTTING
###------------------------------------------

# load term annotation
onto_annot = pd.read_csv(ONTO_ANNOT, sep=";")
latent_terms = onto_annot.Name[onto_annot.depth == 0].tolist()

# load sample annotation
annot = pd.read_csv(SAMPLE_ANNOT)

# create color dict for the annot columns
categs = annot.loc[:,ANNOT_COLS[0]].unique().tolist()
palette = sns.color_palette(cc.glasbey, n_colors=len(categs))
color_dict = dict(zip(categs, palette))

if(len(ANNOT_COLS) > 1):
    categs2 = annot.loc[:,ANNOT_COLS[1]].unique().tolist()
    palette2 = sns.color_palette(cc.glasbey, n_colors=len(categs2))
    color_dict2 = dict(zip(categs2, palette2))

# make scatterplots for all dimensions
if len(latent_terms) % 2 > 0:
    splits = (len(latent_terms) + 1)/2 
else:
    splits = len(latent_terms)/2
dims = np.array_split(range(len(latent_terms)), splits)
if len(latent_terms) % 2 > 0:
    dims[-1] = np.array((dims[-2][1], dims[-1][0]))


for i in dims:
    fig, ax = plt.subplots(1,2, figsize=(20,10))
    sns.scatterplot(x=z[:,i[0]],
            y=z[:,i[1]], 
            hue=annot.loc[:,ANNOT_COLS[0]],
            palette=color_dict,
            legend='full',
            s=8,
            ax=ax[0])
    ax[0].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax[0].set_xlabel(latent_terms[i[0]])
    ax[0].set_ylabel(latent_terms[i[1]])
    if(len(ANNOT_COLS) > 1):
        sns.scatterplot(x=z[:,i[0]],
                y=z[:,i[1]], 
                hue=annot.loc[:,ANNOT_COLS[1]],
                palette=color_dict2,
                legend='full',
                s=8,
                ax=ax[1])
        ax[1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        ax[1].set_xlabel(latent_terms[i[0]])
        ax[1].set_ylabel(latent_terms[i[1]])
    plt.tight_layout()
    run["images/Latent_dims/" + str(i[0])].upload(fig)
    plt.close()








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

# iterate over all hidden layers to plot activations
for d in range(len(mask_list) - 1):

    # get the depth
    act_sub = activation[d].to('cpu').detach().numpy()
    act_sub = np.array(np.split(act_sub, act_sub.shape[1]/PARAMS['neuronnum'], axis=1)).mean(axis=2).T
    #depth = len(mask_list) - 1 - d #encoder
    if act_sub.shape[1] > 1:
        depth = d + 1 #decoder
        depth_terms = onto_annot.Name[onto_annot.depth == depth].tolist()

        # determine dims for plotting
        if len(depth_terms) % 2 > 0:
            splits = (len(depth_terms) + 1)/2 
        else:
            splits = len(depth_terms)/2
        dims = np.array_split(range(len(depth_terms)), splits)
        if len(depth_terms) % 2 > 0:
            dims[-1] = np.array((dims[-2][1], dims[-1][0]))

        # iterate over all terms in depth and plot
        for i in range(len(dims)):
            fig, ax = plt.subplots(1,2, figsize=(20,10))
            sns.scatterplot(x=act_sub[:,dims[i][0]],
                    y=act_sub[:,dims[i][1]], 
                    hue=annot.loc[:,ANNOT_COLS[0]],
                    palette=color_dict,
                    legend='full',
                    s=8,
                    ax = ax[0])
            ax[0].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            ax[0].set_xlabel(depth_terms[dims[i][0]])
            ax[0].set_ylabel(depth_terms[dims[i][1]])
            if len(ANNOT_COLS) > 1:
                sns.scatterplot(x=act_sub[:,dims[i][0]],
                        y=act_sub[:,dims[i][1]], 
                        hue=annot.loc[:,ANNOT_COLS[1]],
                        palette=color_dict2,
                        legend='full',
                        s=8,
                        ax=ax[1])
                ax[1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
                ax[1].set_xlabel(depth_terms[dims[i][0]])
                ax[1].set_ylabel(depth_terms[dims[i][1]])
            plt.tight_layout()
            run["images/depths/" + str(depth) + "/" + str(i)].upload(fig)
            plt.close()
    



run.stop()

