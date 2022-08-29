import numpy as np
import pandas as pd 

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2


project_dir = '/net/data.isilon/ag-cherrmann/ddoncevic/projects/OntoVAE_manuscript/03_manuscript_analysis/04_IFN_response'



# import the gene lists

res = pd.read_csv(project_dir + '/results/PBMC_CD4T_OVAE-119_IFN_signaling_Wilcoxon_results.csv')
predict_genes = res.gene.to_numpy()[0:778]

# get predicted genes that are not in intersection

diff_genes = pd.read_csv(project_dir + '/results/CD4T_IFN-ß_stim_up.txt', header=None)
pred_genes = np.array([g for g in predict_genes if g not in diff_genes.iloc[:,0].tolist()])

# get differential genes for all celltypes
diff_res = pd.read_csv(project_dir + '/results/pbmc_stimulated_diff_genes.csv')

diff_gene = {}
celltypes = diff_res.celltype.unique()


for c in celltypes:
    sub = diff_res[diff_res.celltype == c]
    sub = sub[sub.qval < 0.1]
    sub = sub[sub.stat > 0]
    diff_gene[c] = sub.gene.to_numpy()


# get all celltypes except CD4T

comps = celltypes[celltypes != 'CD4T']


# plot overlap of pred genes with differential genes

diff_gene['CD4T'] = pred_genes

fig, ax = plt.subplots(1,6)
for i in range(len(ax)):
    venn2([set(diff_gene['CD4T']), set(diff_gene[comps[i]])],
            ('CD4T', comps[i]), set_colors=('orange', 'darkgrey'), alpha = 0.8, ax=ax[i])

plt.savefig(project_dir + '/figures/06-01_overlap_pred_but_not_diff_genes_PBMCs.pdf')
plt.clf()