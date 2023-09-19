import pickle
from onto_vae.ontobj import *

project_dir = os.getcwd()


# import the gene lists

res = pd.read_csv(project_dir + '/03_analysis/04_IFN_response/results/Genes_upregulating_IFN_activity_Wilcoxon_results.csv')
predict_genes = res.gene.to_numpy()[0:717]

# get overlap genes

diff_genes = pd.read_csv(project_dir + '/03_analysis/04_IFN_response/results/CD4T_IFN-ß_stim_up.txt', header=None)
ol_genes = np.array([g for g in predict_genes if g in diff_genes.iloc[:,0].tolist()])

# load ontobj
with open(project_dir + '/ontobj/GO_symbol_ontobj.pickle', 'rb') as f:
    go = pickle.load(f) 

# load genes
genes = pd.DataFrame(go.extract_genes())
genes = genes[genes.iloc[:,0].isin(ol_genes)]

# get indices in right order
idx = np.array([genes[genes.iloc[:,0] == o].index.to_numpy()[0] for o in ol_genes])

# load original expression data of control cells
expr = go.extract_dataset('Kang_PBMC_control')

# subset expr data
expr_sub = expr[:,idx]

# convert to dataframe and reshape
expr_data = pd.DataFrame(expr_sub)
expr_data.columns = ol_genes
expr_data = expr_data.melt()
expr_data.columns = ['gene', 'log_tpm']


# make boxplots
fig, ax = plt.subplots(figsize=(4,8))
g = sns.scatterplot(x='log_tpm', y='gene', data=expr_data, color='darkgrey', s=6, rasterized=True)
g = sns.boxplot(x='log_tpm', y='gene', data=expr_data, color='skyblue', showfliers=False)
g.figure.savefig(project_dir + '/03_analysis/04_IFN_response/figures/09-01_expression_of_overlapping_genes.pdf')
plt.clf()

# calculate percentage of cells that express the gene
pct =  np.sum(expr_sub != 0, axis=0)/expr_sub.shape[0]

# make bar plot
pct_data = pd.DataFrame({'gene': ol_genes,
                        'pct_expr': pct})

fig, ax = plt.subplots(figsize=(4,8))
g = sns.barplot(x='pct_expr', y='gene', data=pct_data, color='skyblue')
g.figure.savefig(project_dir + '/03_analysis/04_IFN_response/figures/09-02_pct_expressed_of_overlapping_genes.pdf')
plt.clf()
