
# import python modules
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

# set working directory
working_dir = '/net/data.isilon/ag-cherrmann/ddoncevic/projects/OntoVAE_manuscript'

# load expression data and annot
expr = np.load(working_dir + '/datasets/VEGA_PBMC/data/vega_pbmc_GO_trimmed_expr.npy')
annot = pd.read_csv(working_dir + '/datasets/VEGA_PBMC/data/vega_pbmc_annot.csv')
annot = annot.drop('index', axis=1)
genes = pd.read_csv(working_dir + '/ontologies/GO/genes/GO_symbol_trimmed_genes.txt', header=None)

# get celltypes
cts = annot.cell_type.unique()

# create dict to store cell indices for each cell type
ct_indices = {}
for ct in cts:
    ct_indices[ct] = annot[annot.cell_type == ct].index.to_numpy()


# function to compute wilcoxon test
def compute_stats(idx, celltype, alternative='two-sided'):

    """
    alternative: one of 'two-sided', 'greater', 'less'
    """

    annot_sub = annot.iloc[idx,:]
    ctrl = annot_sub[annot_sub.condition == 'control'].index.to_numpy()
    stim = annot_sub[annot_sub.condition == 'stimulated'].index.to_numpy()

    wilcox = [stats.ranksums(expr[stim,i], expr[ctrl,i]) for i in range(expr.shape[1])]
    stat = np.array([i[0] for i in wilcox])
    pvals = np.array([i[1] for i in wilcox])
    qvals = fdrcorrection(np.array(pvals))

    log2fc = [np.log2(np.median(expr[stim,i])/np.median(expr[ctrl,i])) for i in range(expr.shape[1])]
    res = pd.DataFrame({'celltype': [celltype] * len(pvals),
                        'gene': genes.iloc[:,1].tolist(),
                        'stat': stat,
                        'pval' : pvals,
                        'qval': qvals[1],
                        'log2fc': log2fc})
    return res


# initialize results
res = []

# compute Wilcoxon tests for all cell types separately
for ct in ct_indices:
    res.append(compute_stats(ct_indices[ct], ct))

# concatenate results and save
res = pd.concat(res)
res.to_csv(working_dir + '/03_manuscript_analysis/04_IFN_response/results/pbmc_stimulated_diff_genes.csv', index=False)


# isolate genes upregulated in CD4T cells and save as separate list
res_sig = res[res.qval < 0.1]
res_sig_stim = res_sig[res_sig.log2fc < 0]
res_sig_ct = res_sig_stim[res_sig_stim.celltype == 'CD4T'].sort_values(['qval', 'stat'], ascending=(True, False))
res_sig_ct.gene.to_csv(working_dir + '/03_manuscript_analysis/04_IFN_response/results/CD4T_IFN-ß_stim_up.txt', index=False, header=None)





