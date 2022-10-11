
# import modules and define paths
import pickle
from onto_vae.ontobj import *

project_dir = os.getcwd()


# load pathway activations
act = np.load(project_dir + '/activities/recount3_GTEx_GO_activities_neuronnum3_2.npy')

# import ontobj
with open(project_dir + '/ontobj/GO_ensembl_ontobj.pickle', 'rb') as f:
    go = pickle.load(f) 

# load sample annotation
gtex_annot = pd.read_csv(project_dir + "/datasets/recount3_GTEx/recount3_GTEx_annot.csv")

# helper function to make scatterplots
def make_scatterplot(term1, term2, name):
    go.plot_scatter(sample_annot=gtex_annot,
                        color_by='tissue',
                        act=act,
                        term1=term1,
                        term2=term2)
    plt.savefig(project_dir + '/03_analysis/01_pathway_activities/figures/' + name)
    plt.close()


# example terms
terms = ['digestive system process', 'glutamate receptor signaling pathway', 'aortic valve morphogenesis', 'axon ensheathment']

# make plots
make_scatterplot(terms[0], terms[1], '02-01_GTEx_GO_VAE_scatterplot_example1.pdf')
make_scatterplot(terms[2], terms[3], '02-02_GTEx_GO_VAE_scatterplot_example2.pdf')