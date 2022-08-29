


# import python libraries

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
import seaborn as sns
import colorcet as cc

# define path

working_dir = '/net/data.isilon/ag-cherrmann/ddoncevic/projects/OntoVAE_manuscript'


# dict for the 10 models where link between LGMD and gene was removed

model_dict = {'OVAE-129': 'ACTA1',
              'OVAE-130': 'FHL1',
              'OVAE-131': 'SGCA', 
              'OVAE-132': 'DPM3',
              'OVAE-133': 'TMEM43',
              'OVAE-134': 'LMNA',
              'OVAE-135': 'HNRNPDL',
              'OVAE-136': 'EMD',
              'OVAE-137': 'DAG1',
              'OVAE-138': 'SYNE1'}

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
genes = genes.iloc[:,0].tolist()


# function to rank the genes for one model
def create_res_table(link_model, alternative):

    # load reference model
    z_ref = np.load(working_dir + "/models/" + link_model + "/recount3_GTEx_latent_space_embedding.npy")
    act_ref = np.load(working_dir + "/models/" + link_model + "/recount3_GTEx_HPO_activities.npy")
    act_ref = np.hstack((z_ref, act_ref))
    act_ref = act_ref[muscle_idx,md_ind]

    # initialize results
    wilcox = []

    # iterate over genes and perform diff analysis
    for g in genes:
        print(g)
        act = np.load(working_dir + "/data/03_LGMD/link_removal_knockout/LGMD_" + model_dict[link_model] + "/" + g + "/recount3_GTEx_muscle_" + link_model + "_" + g + "_knockout_LGMD_activities.npy")
        act = act.flatten()
        wilcox.append(stats.wilcoxon(act, act_ref, zero_method='zsplit', alternative=alternative))

    stat = np.array([i[0] for i in wilcox])
    pvals = np.array([i[1] for i in wilcox])
    qvals = fdrcorrection(np.array(pvals))

    res = pd.DataFrame({'link': [model_dict[link_model]] * len(genes),
                        'gene': genes,
                        'stat': stat,
                        'pval': pvals,
                        'qval': qvals[1]})

    return res


# initialize results
res_dn = []
res_up = []

# iterate over models and append results
for model in list(model_dict.keys()):
    res_dn.append(create_res_table(model, 'less'))
    res_up.append(create_res_table(model, 'greater'))

# concatenate 
res_dn_table = pd.concat((res_dn))
res_up_table = pd.concat((res_up))

# save results
res_dn_table.to_csv(working_dir + '/03_manuscript_analysis/03_LGMD/results/removed_link_models_gene_ranking_LGMD_dn.csv')
res_up_table.to_csv(working_dir + '/03_manuscript_analysis/03_LGMD/results/removed_link_models_gene_ranking_LGMD_up.csv')



# load results from full model gene knockouts and filter for their leading edges
ko_res_dn = pd.read_csv(working_dir + '/03_manuscript_analysis/03_LGMD/results/GTEx_muscle_OVAE-120_gene_knockout_LGMD_Wilcoxon_results_dn.csv')
ko_res_up = pd.read_csv(working_dir + '/03_manuscript_analysis/03_LGMD/results/GTEx_muscle_OVAE-120_gene_knockout_LGMD_Wilcoxon_results_up.csv')

# filter for their leading edges
ko_res_dn = ko_res_dn.iloc[0:564,:]
ko_res_up = ko_res_up.iloc[0:424,:]


# determine which LGMD gene is in which leading edge
lgmd_genes = list(model_dict.values())
dn_genes = [l for l in lgmd_genes if l in ko_res_dn.gene.tolist()]
up_genes = [l for l in lgmd_genes if l in ko_res_up.gene.tolist()]


# for each of the 10 genes, find how their ranking changed
rankings = {}

for gene in dn_genes:
    res_sub = res_dn_table[res_dn_table.link == gene].sort_values('qval')
    rankings[gene] = (ko_res_dn.gene.tolist().index(gene) + 1, res_sub.gene.tolist().index(gene) + 1)
for gene in up_genes:
    res_sub = res_up_table[res_up_table.link == gene].sort_values('qval')
    rankings[gene] = (ko_res_up.gene.tolist().index(gene) + 1, res_sub.gene.tolist().index(gene) + 1)



# store in pandas dataframe and save
ranks = pd.DataFrame(rankings).T.reset_index()
ranks.columns = ['gene', 'pre-rank', 'post-rank']
ranks.to_csv(working_dir + '/03_manuscript_analysis/03_LGMD/results/LGMD_link_genes_ranking_pre_and_post_removal.csv', index=False)




# plot the change of the ranks
ranks_up = ranks[ranks.gene.isin(up_genes)]
ranks_dn = ranks[ranks.gene.isin(dn_genes)]

plot_up = pd.melt(ranks_up, id_vars=['gene'], value_vars=['pre-rank', 'post-rank'])
plot_dn = pd.melt(ranks_dn, id_vars=['gene'], value_vars=['pre-rank', 'post-rank'])

ax = sns.scatterplot(x='value', y='gene', hue='variable', data=plot_up)
ax.set_xscale('log')
ax.yaxis.grid(True)
plt.axvline(424, 0, 100)
plt.xlim(0, 3000)
ax.figure.savefig(working_dir + '/03_manuscript_analysis/03_LGMD/figures/05-01_ranking_pre_and_post_link_removal_up_genes.pdf')
plt.clf()

ax = sns.scatterplot(x='value', y='gene', hue='variable', data=plot_dn)
ax.set_xscale('log')
ax.yaxis.grid(True)
plt.axvline(564, 0, 100)
plt.xlim(0, 3000)
ax.figure.savefig(working_dir + '/03_manuscript_analysis/03_LGMD/figures/05-02_ranking_pre_and_post_link_removal_dn_genes.pdf')
plt.clf()