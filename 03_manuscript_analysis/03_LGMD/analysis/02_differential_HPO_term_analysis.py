


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


working_dir = '/net/data.isilon/ag-cherrmann/ddoncevic/projects/OntoVAE_manuscript'



###---------------------> load annotation

# load term annotation
onto_annot = pd.read_csv(working_dir + '/ontologies/HPO/annot/HPO_trimmed_annot.csv', sep=";")
md_ind = onto_annot[onto_annot.Name == 'Limb-girdle muscular dystrophy'].index.to_numpy()

# load gtex annotation
gtex_annot = pd.read_csv(working_dir + "/datasets/recount3_GTEx/data/recount3_GTEx_annot.csv")
gtex_annot_sub = gtex_annot[gtex_annot.tissue == 'Muscle']
muscle_idx = gtex_annot_sub.index.to_numpy()

# load genes
genes = pd.read_csv(working_dir + "/ontologies/HPO/genes/HPO_trimmed_genes.txt", header=None)
genes = genes.iloc[:,0].to_numpy()

# load expr
expr = np.load(working_dir + '/datasets/recount3_GTEx/data/recount3_GTEx_HPO_trimmed_expr.npy')
expr_muscle = expr[muscle_idx,:]

# only keep genes which are not 0 for all cells
good_genes = genes[np.where(np.sum(expr_muscle, axis=0) != 0)[0]]


###-----------------------------------------------------------###
##                     WILCOXON TESTING                        ##
###-----------------------------------------------------------###


# load LGMD activities from trained model
z = np.load(working_dir + '/models/OVAE-120/recount3_GTEx_latent_space_embedding.npy')
act = np.load(working_dir + '/models/OVAE-120/recount3_GTEx_HPO_activities.npy')
act = np.hstack((z,act))
act = act[muscle_idx, md_ind].flatten()


# initialize results
wilcox_up = []
wilcox_dn = []


# iterate over all genes and perform Wilcoxon text knockout vs control
for exp in good_genes:

    # get good cells (cells that are not 0)
    g_idx = np.where(genes == exp)[0]
    good_samp = np.where(expr_muscle[:,g_idx] != 0)[0]  

    ko_act = np.load(working_dir + "/data/03_LGMD/knockout/" + exp + "/recount3_GTEx_muscle_OVAE-120_" + exp + "_knockout_LGMD_activities.npy")
    ko_act = ko_act.flatten()
    wilcox_up.append(stats.wilcoxon(ko_act[good_samp], act[good_samp], alternative='greater')) # genes upregulating the LGMD node
    wilcox_dn.append(stats.wilcoxon(ko_act[good_samp], act[good_samp], alternative='less')) # genes downregulating the LGMD node

w_stat_up = np.array([i[0] for i in wilcox_up])
w_pvals_up = np.array([i[1] for i in wilcox_up])
w_qvals_up = fdrcorrection(np.array(w_pvals_up))
w_stat_dn = np.array([i[0] for i in wilcox_dn])
w_pvals_dn = np.array([i[1] for i in wilcox_dn])
w_qvals_dn = fdrcorrection(np.array(w_pvals_dn))

res_up = pd.DataFrame({'gene': good_genes,
                    'w_stat': w_stat_up,
                    'w_pval': w_pvals_up,
                    'w_qval': w_qvals_up[1]})
res_up = res_up.sort_values('w_pval').reset_index(drop=True)

res_dn = pd.DataFrame({'gene': good_genes,
                    'w_stat': w_stat_dn,
                    'w_pval': w_pvals_dn,
                    'w_qval': w_qvals_dn[1]})
res_dn = res_dn.sort_values('w_pval').reset_index(drop=True)

res_up.to_csv(working_dir + '/03_manuscript_analysis/03_LGMD/results/GTEx_muscle_OVAE-120_gene_knockout_LGMD_Wilcoxon_results_up.csv', index=False)
res_dn.to_csv(working_dir + '/03_manuscript_analysis/03_LGMD/results/GTEx_muscle_OVAE-120_gene_knockout_LGMD_Wilcoxon_results_dn.csv', index=False)

res_up.to_excel(working_dir + '/03_manuscript_analysis/03_LGMD/results/Genes_upregulating_LGMD_activity_Wilcoxon_results.xlsx', index=False)
res_dn.to_excel(working_dir + '/03_manuscript_analysis/03_LGMD/results/Genes_downregulating_LGMD_activity_Wilcoxon_results.xlsx', index=False)





