# This script shows an example of how to prepare a dataset in python for training or running through an OntoVAE model.
# Here, our dataset is provided as h5ad file.
# In the future, this functions for this might be integrated in onto-vae package.
# In order to run this, the neccessary files must have been generated already as shown in 01_preprocess_ontology.py
# The only file we need is a list of genes that could be mapped to the ontology (usually sorted by onto-vae in alphabetical order)
# We will generate two files here:
# 1. 2d numpy array of expression data with samples in rows and features in columns (features will be sorted to map the genes from the ontology)
# 2. a csv-file containing the annotation of the data (not neccessarily needed for training the model but for plotting results)

# import modules
import scanpy
import numpy as np
import pandas as pd

# define directory where dataset is located
data_dir = '/user/'

# load dataset
data = scanpy.read_h5ad(data_dir + '/datasets/scrna_data.h5ad')

# extract elements from scanpy object
expr = data.X.todense()
annot = data.obs
genes = data.var.gene_symbol.reset_index(drop=True)

# convert to pandas dataframe
expr_dat = pd.DataFrame(expr.T)
expr_dat.index = genes
expr_dat.columns = annot.index

# load genes from ontology (this file has been generated in 01_preprocess_ontology.py)
go_genes = pd.read_csv(data_dir + '/ontologies/GO/genes/GO_symbol_trimmed_genes.txt', header=None)
go_genes.index = go_genes.iloc[:,0]

# merge data with GO genes and save
merged_go = go_genes.join(expr_dat).fillna(0).drop(0, axis=1).T
merged_go_np = merged_go.to_numpy()
np.save(data_dir + '/datasets/scrna_data_GO_trimmed_expr.npy', merged_go_np)

# save annotation of data
annot = annot.reset_index(drop=False)
annot.to_csv(data_dir + '/datasets/scrna_data_annot.csv', index=None)

