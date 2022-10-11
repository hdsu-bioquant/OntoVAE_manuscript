
# import modules

import os
import sys
import numpy as np 
import pandas as pd
import json
import pickle
import itertools
import umap
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
import seaborn as sns
import colorcet as cc

project_dir = '/net/data.isilon/ag-cherrmann/ddoncevic/projects/OntoVAE_manuscript'



# import original data

expr = pd.read_csv(project_dir + '/datasets/Kang_PBMC/Kang_PBMC_expr.csv')
expr = expr.to_numpy()
annot = pd.read_csv(project_dir + '/datasets/Kang_PBMC/Kang_PBMC_annot.csv')


# scale the data

scaler = StandardScaler()
expr_scaled = scaler.fit_transform(expr)


# perform PCA

pca = PCA(n_components=30)
expr_pca = pca.fit_transform(expr_scaled)


# perform UMAP

reducer = umap.UMAP()
embedding = reducer.fit_transform(expr_pca)


# create color mapping

categs = annot.loc[:,'celltype'].unique().tolist()
conds = annot.loc[:,'condition'].unique().tolist()
color_dict1 = dict(zip(categs, sns.color_palette(cc.glasbey, n_colors=len(categs))))
color_dict2 = dict(zip(conds, sns.color_palette(cc.glasbey, n_colors=len(conds))))



# create scatter plot

fig, ax = plt.subplots(1,2, figsize=(20,10))
sns.scatterplot(x=embedding[:,0],
                y=embedding[:,1], 
                hue=annot.loc[:,'celltype'],
                palette=color_dict1,
                legend='full',
                s=8,
                ax=ax[0],
                rasterized=True)
ax[0].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
ax[0].set_xlabel('UMAP1')
ax[0].set_ylabel('UMAP2')
sns.scatterplot(x=embedding[:,0],
                y=embedding[:,1], 
                hue=annot.loc[:,'condition'],
                palette=color_dict2,
                legend='full',
                s=8,
                ax=ax[1],
                rasterized=True)
ax[1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
ax[1].set_xlabel('UMAP1')
ax[1].set_ylabel('UMAP2')
for a in ax:
    a.legend([],[], frameon=False)
plt.tight_layout()
plt.savefig(project_dir + '/03_analysis/04_IFN_signaling/figures/01_PBMC_UMAP.pdf')
plt.close()





