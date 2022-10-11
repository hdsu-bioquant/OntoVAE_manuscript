
# import python modules
import pickle
from onto_vae.ontobj import *

project_dir = os.getcwd()

# import ontobj 
with open(project_dir + '/ontobj/GO_symbol_ontobj.pickle', 'rb') as f:
    go = pickle.load(f) 

# match dataset
go.match_dataset(expr_data = project_dir + '/datasets/Kang_PBMC/Kang_PBMC_expr.csv',
                  name='Kang_PBMC')

# extract dataset
expr = go.extract_dataset('Kang_PBMC')

# extract genes
genes = go.extract_genes()

# load sample annotation
annot = pd.read_csv(project_dir + '/datasets/Kang_PBMC/Kang_PBMC_annot.csv')

# get celltypes
cts = annot.celltype.unique()

# create dict to store cell indices for each cell type
ct_indices = {}
for ct in cts:
    ct_indices[ct] = annot[annot.celltype == ct].index.to_numpy()


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

    res = pd.DataFrame({'celltype': [celltype] * len(pvals),
                        'gene': genes,
                        'stat': stat,
                        'pval' : pvals,
                        'qval': qvals[1]})
    return res


# initialize results
res = []

# compute Wilcoxon tests for all cell types separately
for ct in ct_indices:
    res.append(compute_stats(ct_indices[ct], ct))

# concatenate results and save
res = pd.concat(res)
res.to_csv(project_dir + '/03_analysis/04_IFN_response/results/pbmc_stimulated_diff_genes.csv', index=False)


# isolate genes upregulated in CD4T cells and save as separate list
res_sig = res[res.qval < 0.1]
res_sig_stim = res_sig[res_sig.stat > 0] 
res_sig_ct = res_sig_stim[res_sig_stim.celltype == 'CD4T'].sort_values(['qval', 'stat'], ascending=(True, False))
res_sig_ct.gene.to_csv(project_dir + '/03_analysis/04_IFN_response/results/CD4T_IFN-ß_stim_up.txt', index=False, header=None)





