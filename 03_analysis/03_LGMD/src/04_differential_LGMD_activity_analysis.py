import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

project_dir = os.getcwd()



# import LGMD activities
act = np.load(project_dir + '/activities/GTEx_muscle_LGMD_activities.npy')

# import systematic KO LGMD activities
ko_act = pd.read_csv(project_dir + '/03_analysis/03_LGMD/results/systematic_knockout_GTEx_muscle_LGMD_activities.npy')

# initialize results
wilcox_up = []
wilcox_dn = []

# iterate over all genes and perform Wilcoxon text knockout vs control
genes = ko_act.columns.to_numpy()
for gene in genes:

    wilcox_up.append(stats.wilcoxon(ko_act[gene].to_numpy(), act, alternative='greater')) # genes upregulating the LGMD node
    wilcox_dn.append(stats.wilcoxon(ko_act[gene].to_numpy(), act, alternative='less')) # genes downregulating the LGMD node


# function to extract test statistics
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

# extract test stats
res_up = extract_stats(wilcox_up)
res_dn = extract_stats(wilcox_dn)

# write results to file
res_up.to_csv(working_dir + '/03_analysis/03_LGMD/results/Genes_upregulating_LGMD_activity_Wilcoxon_results.csv', index=False)
res_dn.to_csv(working_dir + '/03_analysis/03_LGMD/results/Genes_downregulating_LGMD_activity_Wilcoxon_results.csv', index=False)

res_up.to_excel(working_dir + '/03_analysis/03_LGMD/results/Genes_upregulating_LGMD_activity_Wilcoxon_results.xlsx', index=False)
res_dn.to_excel(working_dir + '/03_analysis/03_LGMD/results/Genes_downregulating_LGMD_activity_Wilcoxon_results.xlsx', index=False)