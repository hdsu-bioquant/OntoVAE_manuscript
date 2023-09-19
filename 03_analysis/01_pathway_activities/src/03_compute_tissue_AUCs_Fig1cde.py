import os
import pickle
from onto_vae.ontobj import *

from sklearn.model_selection import StratifiedKFold
from sklearn.naive_bayes import GaussianNB
from sklearn import metrics

from matplotlib import colors

project_dir = os.getcwd()

# import ontobj
with open(project_dir + '/ontobj/GO_ensembl_ontobj.pickle', 'rb') as f:
    go = pickle.load(f) 

# extract ontology annot
onto_annot = go.extract_annot()
mapping = dict(zip(onto_annot.ID.tolist(), onto_annot.Name.tolist()))

# import pathway activities
act = np.load(project_dir + '/activities/decoder/recount3_GTEx_GO_activities_neuronnum3_2.npy')

# load sample annot
annot = pd.read_csv(project_dir + '/datasets/recount3_GTEx/recount3_GTEx_annot.csv')

# store indices for different tissues
tissues = annot.tissue.unique()
tissue_indices = {}
for tissue in tissues:
    tissue_indices[tissue] = annot[annot.tissue == tissue].index.tolist()



# helper function to calculate auc with 10-fold cross validation
def calc_auc_median(X,y):

    # initialize aucs list
    aucs = []

    # cross-val
    skf = StratifiedKFold(n_splits=10, shuffle=True)

    # iterate over folds
    for train_ind, test_ind in skf.split(X,y,y):
        X_train, X_test, y_train, y_test = X[train_ind], X[test_ind], y[train_ind], y[test_ind]

        # train Naive Bayes
        gnb = GaussianNB()
        gnb.fit(X_train, y_train)

        # make predictions
        y_pred = gnb.predict(X_test)

        # calculate AUC
        fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred)
        auc = metrics.auc(fpr, tpr)

        # append to list
        aucs.append(auc)

    return np.nanmedian(np.array(auc))


# function to calculate all AUCs for a given tissue (1-vs-all classification)
def NB_classify(tissue):
    # get correct 1 vs all labels
    annot['labels'] = annot['tissue'].apply(lambda x: 1 if x == tissue else 0)
    y = annot['labels'].to_numpy()
    # calculate auc for all terms with helper function
    auc = [calc_auc_median(act[:,i].reshape(-1,1), y) for i in range(act.shape[1])]
    auc = pd.DataFrame({tissue: auc})
    return auc



# calculate aucs for all terms and all tissues (1-vs-all)
all_aucs = []

for t in tissues:
    auc = NB_classify(t)
    all_aucs.append(auc)

all_aucs = pd.concat(all_aucs, axis=1)
all_aucs.index = onto_annot.ID

# 10 example terms
example_terms = ['regulation of alpha-beta T cell proliferation',
                 'digestion',
                 'brain development',
                 'aortic valve morphogenesis',
                 'gland development',
                 'muscle system process',
                 'myelination', 
                 'neutral lipid metabolic process',
                 'drug metabolic process',
                 'regulation of vascular permeability']

# create color palette
palette = sns.color_palette(cc.glasbey, n_colors=len(tissues))
color_dict = dict(zip(tissues, palette))

# plot the aucs as dot plot (Fig 1c)

plot_data = all_aucs[all_aucs.index.isin(example_terms)].melt(ignore_index=False)
plot_data = plot_data.reset_index()
plot_data.columns = ['GO term', 'tissue', 'AUC']

fig, ax = plt.subplots(1,1, figsize=(10,10))
ax = sns.stripplot(x="AUC", 
                y="GO term",
                hue="tissue", 
                data=plot_data,
                palette=color_dict,
                s=10,
                alpha=0.9)
ax.yaxis.grid(True)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.tight_layout()
plt.savefig(project_dir + '/03_analysis/01_pathway_activities/figures/03_AUC_stripplot.pdf')
plt.clf()

# make density plots for the ten different terms (Fig 1d)

fig, ax = plt.subplots(10,1, figsize=(5,10))

for i in range(len(example_terms)):

    sub_data = plot_data[plot_data['GO term'] == example_terms[i]]
    sub_tissues = sub_data[sub_data.AUC > 0.5].tissue.to_numpy()

    idx = onto_annot[onto_annot.Name == example_terms[i]].index.to_numpy()[0]
    annot[example_terms[i]] = act[:,idx]
    annot_sub = annot[annot.tissue.isin(sub_tissues)]

    sns.kdeplot(data=annot_sub,
                x=example_terms[i], 
                hue='tissue',
                palette=color_dict,
                fill=True,
                alpha=0.5,
                legend='full',
                ax=ax.flatten()[i])

plt.tight_layout()
plt.savefig(project_dir + '/03_analysis/01_pathway_activities/figures/03_AUC_activity_densities.pdf')
plt.close()




