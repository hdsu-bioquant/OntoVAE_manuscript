
import pickle
from onto_vae.ontobj import *



project_dir = os.getcwd()



# import ontobj
with open(project_dir + '/ontobj/GO_ensembl_ontobj.pickle', 'rb') as f:
    go = pickle.load(f) 


# import activities and reconstructed values
act = np.load(project_dir + '/03_analysis/02_gene_knockout/results/activities/recount3_GTEx_muscle_control_activities.npy')
rec = np.load(project_dir + '/03_analysis/02_gene_knockout/results/reconstructions/recount3_GTEx_muscle_control_reconstructions.npy')


# helper function to perform wilcoxon testing
def wilcox_testing(gene, vartype, option):
    # gene: the ko gene
    # vartype: activities or reconstructions
    # option: terms or genes
    ko = np.load(project_dir + '/03_analysis/02_gene_knockout/results/' + vartype + '/recount3_GTEx_muscle_' + gene + '_knockout_' + vartype + '.npy')

    if vartype == 'activities':
        ctrl = act
    else:
        ctrl = rec    

    res = go.wilcox_test(dataset='recount3_GTEx_muscle',
                         control = ctrl,
                         perturbed = ko,
                         direction='down',
                         option=option
                         )
    return res
    


# get most significant terms
DMD_terms = wilcox_testing('DMD', 'activities', 'terms')
DMD_terms.to_csv(project_dir + '/03_analysis/02_gene_knockout/results/differential_analysis/DMD_ranked_terms.csv', sep=';', index=False)
SFSWAP_terms = wilcox_testing('SFSWAP', 'activities', 'terms')
SFSWAP_terms.to_csv(project_dir + '/03_analysis/02_gene_knockout/results/differential_analysis/SFSWAP_ranked_terms.csv', sep=';', index=False)
COX5A_terms = wilcox_testing('COX5A', 'activities', 'terms')
COX5A_terms.to_csv(project_dir + '/03_analysis/02_gene_knockout/results/differential_analysis/COX5A_ranked_terms.csv', sep=';', index=False)

# get most significant genes
DMD_genes = wilcox_testing('DMD', 'reconstructions', 'genes')
DMD_genes.to_csv(project_dir + '/03_analysis/02_gene_knockout/results/differential_analysis/DMD_ranked_genes.csv', sep=';', index=False)
SFSWAP_genes = wilcox_testing('SFSWAP', 'reconstructions', 'genes')
SFSWAP_genes.to_csv(project_dir + '/03_analysis/02_gene_knockout/results/differential_analysis/SFSWAP_ranked_genes.csv', sep=';', index=False)
COX5A_genes = wilcox_testing('COX5A', 'reconstructions', 'genes')
COX5A_genes.to_csv(project_dir + '/03_analysis/02_gene_knockout/results/differential_analysis/COX5A_ranked_genes.csv', sep=';', index=False)




