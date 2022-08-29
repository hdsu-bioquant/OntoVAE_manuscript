


###---------------------> set up the environment

import sys
import numpy as np 
import pandas as pd
import json
import pickle
import itertools
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

import matplotlib.pyplot as plt
import seaborn as sns
import colorcet as cc


working_dir = '/net/data.isilon/ag-cherrmann/ddoncevic/projects/OntoVAE_manuscript'


###---------------------> load annotation

# load go annotation
go_annot = pd.read_csv(working_dir + '/ontologies/GO/annot/GO_ensembl_trimmed_annot.csv', sep=";")

# load gtex annotation
gtex_annot = pd.read_csv(working_dir + "/datasets/recount3_GTEx/data/recount3_GTEx_annot.csv")
gtex_annot_sub = gtex_annot[gtex_annot.tissue == 'Muscle']
muscle_idx = gtex_annot_sub.index.to_numpy()

# load genes
genes = pd.read_csv(working_dir + '/ontologies/GO/genes/GO_ensembl_trimmed_genes.txt', header=None)
genes = genes.iloc[:,0].to_numpy()



###---------------------> function for Wilcoxon testing

def compute_stats(ctrl, ko, gene, alternative='greater'):

    """
    alternative: one of 'two-sided', 'greater', 'less'
    in this case, greater: genes that are lower in knockout compared to control
    """

    wilcox = [stats.wilcoxon(ctrl[:,i], ko[:,i], zero_method='zsplit', alternative=alternative) for i in range(ctrl.shape[1])]
    stat = np.array([i[0] for i in wilcox])
    pvals = np.array([i[1] for i in wilcox])
    qvals = fdrcorrection(np.array(pvals))

    res = pd.DataFrame({'ko_gene': gene,
                        'gene': genes,
                        'stat': stat,
                        'pval' : pvals,
                        'log10p': -np.log10(pvals),
                        'qval': qvals[1]})
    return(res)



###---------------------> compute Wilcoxon test between knockout and control

res = []

rec = np.load(working_dir + "/models/OVAE-113/recount3_GTEx_reconstruction.npy")
rec = rec[muscle_idx,:]


for gene in  ['DMD', 'SFSWAP', 'COX5A']:
    ko_rec = np.load(working_dir + "/data/02_gene_knockout/" + gene + "/recount3_GTEx_muscle_OVAE-113_" + gene + "_knockout_reconstruction.npy")
    res.append(compute_stats(rec, ko_rec, gene))

res = pd.concat(res)

res.to_csv(working_dir + '/03_manuscript_analysis/02_gene_knockout/results/GTEx_muscle_OVAE-113_gene_knockout_reconstruction_Wilcoxon_results.csv', sep=";", index=False)
res.to_excel(working_dir + '/03_manuscript_analysis/02_gene_knockout/results/GTEx_muscle_OVAE-113_gene_knockout_reconstruction_Wilcoxon_results.xlsx', index=False)


