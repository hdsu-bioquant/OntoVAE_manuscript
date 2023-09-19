import os
import numpy as np
import pandas as pd
import pickle
import umap

from onto_vae.ontobj import *
from onto_vae.vae_model import *

project_dir = os.getcwd()

# import ontobj
with open(project_dir + '/ontobj/GO_symbol_ontobj.pickle', 'rb') as f:
    go = pickle.load(f) 


# initialize OntoVAE model
go_model = OntoVAE(ontobj=go,
                dataset='Kang_PBMC_control',
                top_thresh=1000,
                bottom_thresh=30)
go_model.to(go_model.device)


# load best model
checkpoint = torch.load(project_dir + '/models/PBMC_GO/best_model.pt', 
                        map_location=torch.device(go_model.device))
go_model.load_state_dict(checkpoint['model_state_dict'])

# get the ontology genes
onto_genes = go.extract_genes()

# load ranked gene list for IFN stimulation (leading edge is at 717)
IFN_genelist = pd.read_csv(project_dir + '/03_analysis/04_IFN_response/results/Genes_upregulating_IFN_activity_Wilcoxon_results.csv')

# get gene indices
IFN_gene_inds = [onto_genes.index(g) for g in IFN_genelist.gene.to_numpy()[:717]]
ctrl_inds =  [onto_genes.index(g) for g in IFN_genelist.gene.to_numpy()[::-1][:717]]

# create cutoffs at 10,20,... %
top_genelists = [IFN_genelist.gene.to_numpy()[:int(np.round(717*p))] for p in np.arange(0.1, 1.1, 0.1)]
top_ind = [IFN_gene_inds[:int(np.round(717*p))] for p in np.arange(0.1, 1.1, 0.1)]

# load sample annotation
sample_annot =  pd.read_csv(project_dir + '/datasets/Kang_PBMC/Kang_PBMC_annot.csv')
cd4t_ctrl_idx = sample_annot[(sample_annot.celltype == 'CD4T') & (sample_annot.condition == 'control')].index.to_numpy()
cd4t_stim_idx = sample_annot[(sample_annot.celltype == 'CD4T') & (sample_annot.condition == 'stimulated')].index.to_numpy()

# extract dataset
data = go.extract_dataset('Kang_PBMC',
                            top_thresh=1000,
                            bottom_thresh=30)


# add subsetted datasets to Ontoobj
go.add_dataset(dataset = data[cd4t_ctrl_idx,:],
               name = 'CD4T_PBMC_control',
               top_thresh=1000,
               bottom_thresh=30)
go.add_dataset(dataset = data[cd4t_stim_idx,:],
               name = 'CD4T_PBMC_stimulated',
               top_thresh=1000,
               bottom_thresh=30)


# get activities of control cells
ctrl_act = go_model.get_pathway_activities(go, 'CD4T_PBMC_control')
np.save(project_dir + '/03_analysis/04_IFN_response/results/CD4T_ctrl_activities.npy', ctrl_act)

# get activities of stim cells
stim_act = go_model.get_pathway_activities(go, 'CD4T_PBMC_stimulated')
np.save(project_dir + '/03_analysis/04_IFN_response/results/CD4T_IFN_stim_activities.npy', stim_act)


# function to stimulate cells in silico
def stim_ctrl_data(gene_inds, operation):
    ctrl_data = go.extract_dataset('CD4T_PBMC_control')
    for i in gene_inds:
        non_zeros = np.where(ctrl_data[:,i] != 0)[0]
        ctrl_data[non_zeros,i] = ctrl_data[non_zeros,i] + operation
    return ctrl_data
 

# stimulate only x % perturbed, but with increasing numbers
silicox_act = []

rand_ind = np.random.choice(np.array(IFN_gene_inds), size=50, replace=False)

for op in np.arange(2, 10, 2):
    stim_data = stim_ctrl_data(rand_ind,operation=op)
    go.add_dataset(dataset = stim_data,
                    name = 'CD4T_in_silico')
    act = go_model.get_pathway_activities(go, 'CD4T_in_silico')
    silicox_act.append(act)   


# create annotation data for plot
plot_annot = pd.DataFrame({'group': ['ground_truth'] * len(cd4t_stim_idx) + 
                                    ['control'] * len(cd4t_ctrl_idx) + 
                                    list(np.repeat([str(p) + ' added' for p in np.arange(2,10,2)], len(cd4t_ctrl_idx)))})

# stack data together
act_data = np.vstack([stim_act, ctrl_act] + silicox_act)

# run UMAP on original data and project in silico stimulated ones
reducer = umap.UMAP(random_state=62)
embedding1 = reducer.fit_transform(act_data[:len(cd4t_stim_idx) + len(cd4t_ctrl_idx),:])
embedding2 = reducer.transform(act_data[len(cd4t_stim_idx) + len(cd4t_ctrl_idx):,:])
embedding = np.vstack([embedding1, embedding2])


# make scatterplots with density plots

for i in range(len(stims)):

    annot_sub = plot_annot[plot_annot.group.isin(['ground_truth', 'control', str(stims[i]) + ' added'])]    

    covar_categs = annot_sub['group'].unique().tolist()
    palette = sns.color_palette(cc.glasbey, n_colors=len(covar_categs))
    color_dict = dict(zip(covar_categs, palette))

    g = sns.JointGrid(x=embedding[annot_sub.index.to_numpy(),0].flatten(),
                    y=embedding[annot_sub.index.to_numpy(),1].flatten(), 
                    hue=annot_sub['group'],
                    palette=color_dict,
                    xlim=(-5,2))

    g.plot_joint(sns.scatterplot, s=15, rasterized=True)
    g.plot_marginals(sns.kdeplot, fill=True,
                    alpha=0.5,)
    plt.savefig(project_dir + '/03_analysis/04_IFN_response/figures/06_leading_edge_UMAP/UMAP_density_add' + str(stims[i]) + '.pdf')
    plt.clf()



# make density plots only
stims = np.arange(2,10,2)

fig, ax = plt.subplots(4,1, figsize=(5,20))

for i in range(len(stims)):

    annot_sub = plot_annot[plot_annot.group.isin(['ground_truth', 'control', str(stims[i]) + ' added'])]

    #create color dict
    covar_categs = annot_sub['group'].unique().tolist()
    palette = sns.color_palette(cc.glasbey, n_colors=len(covar_categs))
    color_dict = dict(zip(covar_categs, palette))

    # make scatter plot
    sns.kdeplot(y=embedding[annot_sub.index.to_numpy(),1].flatten(), 
                    hue=annot_sub['group'],
                    palette=color_dict,
                    fill=True,
                    alpha=0.5,
                    legend='full',
                   ax=ax.flatten()[i])

plt.tight_layout()
plt.savefig(project_dir + '/03_analysis/04_IFN_response/figures/06-02_rand50_leading_edge_genes_UMAP_density.pdf')
plt.close()