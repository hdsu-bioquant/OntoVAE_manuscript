# import python libraries

import os

import numpy as np 
import pandas as pd
import json
import pickle
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

import matplotlib.pyplot as plt
import seaborn as sns
import colorcet as cc

# define path
working_dir = os.getcwd()


# helper function to extract test stats
def extract_stats(wilcox):
    w_stat = np.array([i[0] for i in wilcox])
    w_pval =  np.array([i[1] for i in wilcox])
    w_qval = fdrcorrection(np.array(w_pval))
    res = pd.DataFrame({'gene':genes,
                            'w_stat': w_stat,
                            'w_pval': w_pval,
                            'w_qval': w_qval[1]})
    res = res.sort_values('w_pval').reset_index(drop=True)
    return res


# function to rank genes for one model (basically same as script #04)
def create_res_table(act_path, ko_act_path, alternative):

    # load activities
    act = np.load(act_path)
    ko_act = pd.read_csv(ko_act_path)

    # initialize results
    wilcox = []

    # iterate over genes and perform Wilcoxon test ko vs ctrl
    genes = ko_act.columns.to_numpy()
    for gene in genes:
        wilcox.append(stats.wilcoxon(ko_act[gene].to_numpy(), act, alternative=alternative)) 
   
    # extract test stats
    res = extract_stats(wilcox)
    return res


# load the LGMD link removal genes
lm_genes = pd.read_csv(project_dir + '/03_analysis/03_LGMD/results/LGMD_genes_in_leading_edge.csv')


# initialize results
res_dn = []
res_up = []

# iterate over models and append results
for gene in lm_genes.iloc[:,0].tolist():
    act = project_dir + '/03_analysis/03_LGMD/results/link_removal/GTEx_muscle_LGMD_activities_' + gene + '_link_removed.npy'
    ko_act = project_dir + '/03_analysis/03_LGMD/results/link_removal/systemtic_knockout_GTEx_muscle_LGMD_activities_' + gene + '_link_removed.npy'
    res_dn.append(create_res_table(act, ko_act, 'less')) # genes downregulating LGMD node
    res_up.append(create_res_table(act, ko_act, 'greater')) # genes upregulating LGMD node

# concatenate 
res_dn_table = pd.concat((res_dn))
res_up_table = pd.concat((res_up))

# save results
res_dn_table.to_csv(working_dir + '/03_analysis/03_LGMD/results/link_removal/removed_link_models_gene_ranking_LGMD_dn.csv')
res_up_table.to_csv(working_dir + '/03_analysis/03_LGMD/results/link_removal/removed_link_models_gene_ranking_LGMD_up.csv')



# load results from full model gene knockouts and filter for their leading edges
ko_res_dn = pd.read_csv(working_dir + '/03_analysis/03_LGMD/results/Genes_downregulating_LGMD_activity_Wilcoxon_results.csv')
ko_res_up = pd.read_csv(working_dir + '/03_analysis/03_LGMD/results/Genes_upregulating_LGMD_activity_Wilcoxon_results.csv')

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
ranks.to_csv(working_dir + '/03_analysis/03_LGMD/results/link_removal/LGMD_link_genes_ranking_pre_and_post_removal.csv', index=False)




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
ax.figure.savefig(working_dir + '/03_analysis/03_LGMD/figures/11-01_ranking_pre_and_post_link_removal_up_genes_Fig3e.pdf')
plt.clf()

ax = sns.scatterplot(x='value', y='gene', hue='variable', data=plot_dn)
ax.set_xscale('log')
ax.yaxis.grid(True)
plt.axvline(564, 0, 100)
plt.xlim(0, 3000)
ax.figure.savefig(working_dir + '/03_analysis/03_LGMD/figures/11-02_ranking_pre_and_post_link_removal_dn_genes_Fig3d.pdf')
plt.clf()