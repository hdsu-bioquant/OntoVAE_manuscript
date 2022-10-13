import os
import pickle
from onto_vae.ontobj import *

project_dir = os.getcwd()

# import ontobj
with open(project_dir + '/ontobj/GO_ensembl_ontobj.pickle', 'rb') as f:
    go = pickle.load(f) 

# extract ontology annot
onto_annot = go.extract_annot()

# import pathway activities
act = np.load(project_dir + '/activities/recount3_GTEx_GO_activities_neuronnum3_2.npy')

# load sample annot
annot = pd.read_csv(project_dir + '/datasets/recount3_GTEx/recount3_GTEx_annot.csv')

# store indices for different tissues
tissues = annot.tissue.unique()
tissue_indices = {}
for tissue in tissues:
    tissue_indices[tissue] = annot[annot.tissue == tissue].index.tolist()


# function to calculate all pairwise tissue comparisons for a term and summarize results in table
def compute_stats(term):

    act_sub = act[:,term]

    hits = []
    med_stat = []

    for obj in tissues:
        ref = [t for t in tissues if t != obj]
        wilcox = [stats.ranksums(act_sub[tissue_indices[obj]], act_sub[tissue_indices[r]]) for r in ref]
        stat = np.array([i[0] for i in wilcox])
        pval = np.array([i[1] for i in wilcox])
        hits.append(np.sum(np.logical_and(stat > 0, pval < 0.05)))
        med_stat.append(np.median(stat[np.logical_and(stat > 0, pval < 0.05)]))
    
    table = pd.DataFrame({'ind': [term] * len(tissues),
                          'id': [onto_annot.iloc[term,:].ID] * len(tissues),
                          'term': [onto_annot.iloc[term,:].Name] * len(tissues),
                        'genes': [onto_annot.iloc[term,:].desc_genes] * len(tissues),
                        'tissue': tissues,
                        'hits': hits,
                        'med_stat': med_stat})

    table = table.sort_values(['hits', 'med_stat'], ascending=False)
    table['rank'] = range(1,len(tissues) + 1)
    return table


# iterate over all terms and append results
results = []

for term in range(act.shape[1]):
    results.append(compute_stats(term))

results = pd.concat(results)

# save the results table
results.to_csv(project_dir + '/03_analysis/01_pathway_activities/results/GTEx_tissues_top_GO_terms.csv', index=False, sep=';')

# as sheet for supplementary table
with pd.ExcelWriter(project_dir + '/03_analysis/01_pathway_activities/results/GTEx_tissues_top_GO_terms.xlsx') as writer:  
    for tissue in tissues:
        res_sub = results[results.tissue == tissue]
        res_sub = res_sub.sort_values(['rank', 'hits', 'med_stat'], ascending = (True,False,False))
        res_sub = res_sub.drop(['ind'], axis=1)
        res_sub.to_excel(writer, sheet_name=tissue, index=False)


