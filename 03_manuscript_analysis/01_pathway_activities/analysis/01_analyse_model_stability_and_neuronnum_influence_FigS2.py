# We want to see if the model generates reproducible results when run multiple times with the same parameters.
# We also want to investigate this dependent on the number of neurons per term used.



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import json

project_dir = '/net/data.isilon/ag-cherrmann/ddoncevic/projects/OntoVAE_manuscript'


# define model dict: keys -> neuronnum, values -> model runs

model_dict = {'1': ['OVAE-109', 'OVAE-110'],
              '2': ['OVAE-111', 'OVAE-112'],
              '3': ['OVAE-113', 'OVAE-114']}

runs = [item for sublist in model_dict.values() for item in sublist]

# import latent space embedding and decoder activities for each model
model_act = {}

for run in runs:
    lat = np.load(project_dir + '/models/' + run + '/recount3_GTEx_latent_space_embedding.npy')
    act = np.load(project_dir + '/models/' + run + '/recount3_GTEx_GO_activities.npy')
    model_act[run] = np.hstack((lat,act))

term_no = model_act['OVAE-109'].shape[1]  #3245


# calculate correlations between same neuronnums
pair_corrs = {}

for nn in ['1', '2', '3']:
    pair_corrs[nn] = [np.corrcoef(model_act[model_dict[nn][0]][:,i], model_act[model_dict[nn][1]][:,i])[0,1] for i in range(term_no)]


# calculate cross-correlations between different neuronnums
combos = [('1', '2'), ('1', '3'), ('2', '3')]

cross_corrs = {}

for c in combos:
    cross_corrs[c] = [np.corrcoef(model_act[model_dict[c[0]][0]][:,i], model_act[model_dict[c[1]][0]][:,i])[0,1] for i in range(term_no)]


# make a boxplot of the correlations - Fig S2a
pair_corrs.update(cross_corrs)
flierprops = dict(marker='o', 
                  markersize=3,
                  linestyle='none')

fig, ax = plt.subplots(figsize=(8,4))
g = sns.boxplot(data=list(pair_corrs.values()), color='skyblue', flierprops=flierprops)
g.set_xticklabels(list(pair_corrs.keys()))
g.figure.savefig(project_dir + '/figures/01-01_activity_correlations_based_on_neuronnum.pdf')
plt.clf()


# compute all pairwise correlations
def compute_pathway_corrs(act):
    corrs = {}
    for i in range(act.shape[1]):
        for j in range(i + 1, act.shape[1]):
            corrs[(i,j)] = np.corrcoef(act[:,i], act[:,j])[0,1]
    return corrs

corrs_neuron1 = compute_pathway_corrs(model_act['OVAE-109'])
corrs_neuron2 = compute_pathway_corrs(model_act['OVAE-111'])
corrs_neuron3 = compute_pathway_corrs(model_act['OVAE-113'])

# store results in dict and save them
neuron_corrs = {'1': list(corrs_neuron1.values()),
                '2': list(corrs_neuron2.values()),
                '3': list(corrs_neuron3.values())}

with open(project_dir + "/results/GO_activity_correlations_based_on_neuronnum.json", 'w') as jfile:
    json.dump(neuron_corrs, jfile, sort_keys=True, indent=4)



# plot correlations for upper thresholds

threshs = [0.9, 0.95, 0.99]

counts = [[np.sum(np.array(neuron_corrs[k]) > t) for k in neuron_corrs.keys()] for t in threshs]
counts = pd.DataFrame(counts)

counts.columns = [str(t) for t in threshs]
counts['neuronnum'] = neuron_corrs.keys()
counts = counts.melt(id_vars=['neuronnum'])
counts.columns = ['neuronnum', 'threshold', 'count']


# Fig S2b

g = sns.barplot(x='threshold',y='count', hue='neuronnum', data=counts)
g.figure.savefig(project_dir + '/figures/01-02_Onto_decoder_GO_activity_correlations_based_on_neuronnum_higher_than_thresh.pdf')