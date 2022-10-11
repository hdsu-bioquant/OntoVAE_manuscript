
# set up the environment

import os
import numpy as np 
import pandas as pd
import json


from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import seaborn as sns
import colorcet as cc

from onto_vae.ontobj import *

working_dir = os.getcwd()


# import results of paired Wilcoxon testing 
res_up = pd.read_csv(working_dir + '/03_analysis/03_LGMD/results/Genes_upregulating_LGMD_activity_Wilcoxon_results.csv', index_col=0)
res_dn = pd.read_csv(working_dir + '/03_analysis/03_LGMD/results/Genes_downregulating_LGMD_activity_Wilcoxon_results.csv',index_col=0)

# import DGEA results and filter to downregulated genes
diff_res = pd.read_csv(working_dir + '/03_analysis/03_LGMD/results/LGMD_DGEA_results.csv', index_col=0)
diff_res = diff_res[diff_res.index.isin(res_up.gene.tolist())]
diff_res_sig = diff_res[diff_res.padj < 0.1]
diff_genes_dn = diff_res_sig[diff_res_sig.log2FoldChange < 0].index.to_numpy()

# load LGMD genes in leading edge
lgmd_genes = pd.read_csv(working_dir + '/03_analysis/03_LGMD/results/LGMD_genes_in_leading_edge.csv', header=None)
lgmd_genes = lgmd_genes.iloc[:,0].to_numpy()


# load LGMD annotated genes
# import ontobj
with open(project_dir + '/ontobj/HPO_ontobj.pickle', 'rb') as f:
    hpo = pickle.load(f) 
onto_genes = hpo.desc_genes['1000_10']

lgmd = np.array(onto_genes['HP:0006785'])  # ID for Limb-girdle muscular dystrophy


# make Fig3c (Hockeystick plot)

res_plot = pd.concat((res_up.iloc[0:424,:], res_dn.iloc[0:564]))   # 424 and 564 are the ranks from GSEA
res_plot['log10p'] = -np.log10(res_plot['w_pval'])
res_plot['loglog10p'] = np.log2(res_plot['log10p'])
res_plot = res_plot.sort_values('loglog10p', ascending=False)
res_plot['rank'] = np.arange(1,res_plot.shape[0]+1)

g = sns.lineplot(data=res_plot,
                 x=res_plot['rank'],
                 y=res_plot['loglog10p'],
                 color='black')
g.set(ylabel='log2(-log10(p-value))')


g = sns.scatterplot(data=res_plot[res_plot.gene.isin(lgmd_genes)],
                    x='rank',
                    y='loglog10p',
                    color='darkred')
plt.legend([],[], frameon=False)

for i, g in enumerate(res_plot.gene.to_list()):
    if g in lgmd_genes:
        plt.text(i + 1, 
                res_plot['loglog10p'].iloc[i], 
                g, 
                color='darkgrey') 

plt.savefig(working_dir + '/03_analysis/03_LGMD/figures/06-01_ranking_of_HPO_knockout_predicted_genes_Fig3c.pdf')
plt.clf()


## make Venn diagram with lgmd genes

venn2([set(lgmd), set(res_plot.gene.tolist())], ('lgmd genes', 'predicted genes'))
plt.savefig(working_dir + '/03_analysis/03_LGMD/figures/06-02_overlap_lgmd_genes_leading_edge.pdf')
plt.clf()


## make Venn diagram with differential downregulated genes

venn2([set(diff_genes_dn), set(res_plot.gene.tolist())], ('diff genes dn', 'predicted genes'))
plt.savefig(working_dir + '/03_analysis/03_LGMD/figures/06-03_overlap_diff_genes_dn_leading_edge.pdf')
plt.clf()
