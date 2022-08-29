


###---------------------> set up the environment

import os
import sys
import numpy as np 
import pandas as pd
import json
import pickle
import itertools
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import seaborn as sns
import colorcet as cc


working_dir = '/net/data.isilon/ag-cherrmann/ddoncevic/projects/OntoVAE_manuscript'


###---------------------> load annotation

# load go annotation
onto_annot = pd.read_csv(working_dir + '/ontologies/GO/annot/GO_ensembl_trimmed_annot.csv', sep=";")

# load gtex annotation
gtex_annot = pd.read_csv(working_dir + '/datasets/recount3_GTEx/data/recount3_GTEx_annot.csv')
gtex_annot_sub = gtex_annot[gtex_annot.tissue == 'Muscle']
muscle_idx = gtex_annot_sub.index.to_numpy()

# ko genes
ko_genes = ['DMD', 'SFSWAP', 'COX5A']

###-----------------------------------------------------------###
##                     WILCOXON TESTING                        ##
###-----------------------------------------------------------###

###---------------------> function for Wilcoxon testing

def compute_stats(ctrl, ko, gene, alternative='greater'):

    """
    alternative: one of 'two-sided', 'greater', 'less'
    """

    wilcox = [stats.wilcoxon(ctrl[:,i], ko[:,i], zero_method='zsplit', alternative=alternative) for i in range(ctrl.shape[1])]
    stat = np.array([i[0] for i in wilcox])
    pvals = np.array([i[1] for i in wilcox])
    qvals = fdrcorrection(np.array(pvals))

    #log2fc = [np.log2(np.median(ctrl[:,i])/np.median(ko[:,i])) for i in range(ctrl.shape[1])]
    res = pd.DataFrame({'gene': gene,
                        'id': onto_annot.ID.tolist(),
                        'term': onto_annot.Name.tolist(),
                        'depth': onto_annot.depth.tolist(),
                        'stat': stat,
                        'pval' : pvals,
                        'log10p': -np.log10(pvals),
                        'qval': qvals[1]})#,
                        #'log2fc': log2fc})
    return(res)





###---------------------> compute Wilcoxon test between knockout and control

res = []

z = np.load(working_dir + "/models/OVAE-113/recount3_GTEx_latent_space_embedding.npy")
act = np.load(working_dir + "/models/OVAE-113/recount3_GTEx_GO_activities.npy")
act = np.hstack((z,act))

act = act[muscle_idx,:]
#act[act < 0] = 0

for g in ko_genes:
    ko_z = np.load(working_dir + "/data/02_gene_knockout/" + g + "/recount3_GTEx_muscle_OVAE-113_" + g + "_knockout_latent_space_embedding.npy")
    ko_act = np.load(working_dir + "/data/02_gene_knockout/" + g + "/recount3_GTEx_muscle_OVAE-113_" + g + "_knockout_GO_activities.npy")
    ko_act = np.hstack((ko_z, ko_act))
    #ko_act[ko_act < 0] = 0
    w_stats = compute_stats(act, ko_act, g)
    res.append(w_stats)

res = pd.concat(res)

res.to_csv(working_dir + '/03_manuscript_analysis/02_gene_knockout/results/GTEx_muscle_OVAE-113_gene_knockout_GO_decoder_Wilcoxon_results.csv')









