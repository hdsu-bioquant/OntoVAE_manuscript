# OVAE-8: GO model trained on VEGA PBMC

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

from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import seaborn as sns
import colorcet as cc


project_dir = '/net/data.isilon/ag-cherrmann/ddoncevic/projects/OntoVAE_manuscript'



###---------------------> load annotation

# load go annotation
go_annot = pd.read_csv(project_dir + '/ontologies/GO/annot/GO_symbol_trimmed_annot.csv', sep=";")
inf_ind = go_annot[go_annot.Name == 'type I interferon signaling pathway'].index.to_numpy()

# load sample annotation and subset to unstimulated CD4T
sample_annot = pd.read_csv(project_dir + "/datasets/VEGA_PBMC/data/vega_pbmc_annot.csv")
sample_annot_sub = sample_annot[sample_annot.condition == 'control'].reset_index(drop=True)

sample_annot_sub = sample_annot_sub[sample_annot_sub.celltype == 'CD4T']
cell_idx = sample_annot_sub.index.to_numpy()

# load genes
genes = pd.read_csv(project_dir + '/ontologies/GO/genes/GO_symbol_trimmed_genes.txt', header=None)
genes = genes.iloc[:,0].to_numpy()

# load original expression data
expr = np.load(project_dir + '/datasets/VEGA_PBMC/data/vega_pbmc_control_GO_trimmed_expr.npy')

# filter for cd4t cells
cd4t_expr = expr[cell_idx,:]

# only keep genes which are not 0 for all cells
good_genes = genes[np.where(np.sum(cd4t_expr, axis=0) != 0)[0]]



###---------------------> compute Wilcoxon test between stimulation and control

z = np.load(project_dir + "/models/OVAE-119/VEGA_PBMC_control_latent_space_embedding.npy")
act = np.load(project_dir + "/models/OVAE-119/VEGA_PBMC_control_GO_activities.npy")
act = np.hstack((z,act))
act = act[cell_idx,inf_ind]

wilcox = []

for g in good_genes:

    # get good cells (cells that are not 0)
    g_idx = np.where(genes == g)[0]
    good_cells = np.where(cd4t_expr[:,g_idx] != 0)[0]  

    stim_act = np.load(project_dir + "/data/04_IFN_response/max_stimulation/" + g + "/vega_pbmc_control_OVAE-119_" + g + "_stimulation_interferon_signaling_activity.npy")
    stim_act = stim_act.flatten()
    stim_act = stim_act[cell_idx]
    w_stats = stats.wilcoxon(stim_act[good_cells], act[good_cells], alternative='greater') # higher activity at IFN node
    wilcox.append(w_stats)


w_stat = np.array([i[0] for i in wilcox])
w_pvals = np.array([i[1] for i in wilcox])
w_qvals = fdrcorrection(np.array(w_pvals))


res = pd.DataFrame({'gene': good_genes,
                    'w_stat': w_stat,
                    'w_pval': w_pvals,
                    'w_qval': w_qvals[1]})
res = res.sort_values('w_pval').reset_index(drop=True)


res.to_csv(project_dir + '/03_manuscript_analysis/04_IFN_response/results/PBMC_CD4T_OVAE-119_IFN_signaling_Wilcoxon_results.csv', index=False)
res.to_excel(project_dir + '/03_manuscript_analysis/04_IFN_response/results/Genes_upregulating_IFN_activity_Wilcoxon_results.xlsx', index=False)







