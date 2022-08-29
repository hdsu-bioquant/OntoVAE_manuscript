
# import modules and define paths

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import json

import matplotlib.pyplot as plt
import seaborn as sns
import colorcet as cc

working_dir = '/net/data.isilon/ag-cherrmann/ddoncevic/projects/OntoVAE_manuscript'


# load model data

lat = np.load(working_dir + '/models/OVAE-113/recount3_GTEx_latent_space_embedding.npy')
act = np.load(working_dir + '/models/OVAE-113/recount3_GTEx_GO_activities.npy')
act = np.hstack((lat,act))


# load ontology and sample annot

go_annot = pd.read_csv(working_dir + '/ontologies/GO/annot/GO_ensembl_trimmed_annot.csv', sep=";")
gtex_annot = pd.read_csv(working_dir + "/recount3_GTEx/data/recount3_GTEx_annot.csv")

categs = gtex_annot.loc[:,'tissue'].unique().tolist()
palette = sns.color_palette(cc.glasbey, n_colors=len(categs))
color_dict = dict(zip(categs, palette))


# get indices of terms of interest 

terms = ['digestive system process', 'glutamate receptor signaling pathway', 'aortic valve morphogenesis', 'axon ensheathment']



# function to make scatterplot

def scatter_plot(term1, term2, name):
        ind1 = go_annot[go_annot.Name == term1].index.to_numpy()
        ind2 = go_annot[go_annot.Name == term2].index.to_numpy()

        fig, ax = plt.subplots(figsize=(10,7))
        sns.scatterplot(x=act[:,ind1].flatten(),
                        y=act[:,ind2].flatten(),
                        hue=gtex_annot.loc[:,'tissue'],
                        palette=color_dict,
                        legend='full',
                        s=8,
                        rasterized=True)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.xlabel(term1)
        plt.ylabel(term2)
        plt.tight_layout()
        plt.savefig(working_dir + '/03_manuscript_analysis/01_pathway_activities/figures/' + name)
        plt.close()


# draw the scatterplots

scatter_plot(terms[0], terms[1], '02-01_GTEx_GO_VAE_scatterplot_example1.pdf')
scatter_plot(terms[2], terms[3], '02-02_GTEx_GO_VAE_scatterplot_example2.pdf')