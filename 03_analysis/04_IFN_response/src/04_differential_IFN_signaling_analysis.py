import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection


project_dir = os.getcwd()


# import IFN activities
act = np.load(project_dir + '/activities/PBMC_control_CD4T_IFN_activities.npy')

# import systematic stimulation LGMD activities
stim_act = pd.read_csv(project_dir + '/03_analysis/04_IFN_response/results/systematic_stimulation_control_CD4T_IFN_activities.csv')

# import ontobj 
with open(project_dir + '/ontobj/GO_symbol_ontobj.pickle', 'rb') as f:
    go = pickle.load(f) 

# load Kang PBMC annotation and get CD4T indices
sample_annot = pd.read_csv(project_dir + '/datasets/Kang_PBMC/Kang_PBMC_control_annot.csv')
cd4t_idx = sample_annot[sample_annot.celltype == 'CD4T'].index.to_numpy()

# extract dataset and subset to CD4T
data = go.extract_dataset('Kang_PBMC_control',
                        top_thresh=1000,
                        bottom_thresh=30)
cd4t_data = data[cd4t_idx,:]


# extract ontology genes
onto_genes = go.extract_genes(top_thresh=1000,
                                bottom_thresh=30)

# initialize results
wilcox = []
n_cells = []

# iterate over all genes and perform Wilcoxon text knockout vs control
for gene in good_genes:
    gene_idx = onto_genes.index(gene)
    good_cells = np.where(cd4t_data[:,gene_idx] != 0)[0]
    n_cells.append(len(good_cells))   
    wilcox.append(stats.wilcoxon(stim_act[gene].iloc[good_cells].to_numpy(), act[good_cells],  alternative='greater')) # genes upregulating the IFN node
    

# function to extract test statistics
def extract_stats(wilcox):
    w_stat = np.array([i[0] for i in wilcox])
    w_pval =  np.array([i[1] for i in wilcox])
    w_qval = fdrcorrection(np.array(w_pval))
    res = pd.DataFrame({'gene': good_genes,
                        'w_stat': w_stat,
                        'n': n_cells,
                        'w_pval': w_pval,
                        'w_qval': w_qval[1]})
    res = res.sort_values('w_pval').reset_index(drop=True)
    return res


# extract test stats
res = extract_stats(wilcox)

# write results to file
res.to_csv(project_dir + '/03_analysis/04_IFN_response/results/Genes_upregulating_IFN_activity_Wilcoxon_results.csv', index=False)
res.to_excel(project_dir + '/03_analysis/04_IFN_response/results/Genes_upregulating_IFN_activity_Wilcoxon_results.xlsx', index=False)







