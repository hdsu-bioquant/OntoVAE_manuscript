# This runs the script train_model.py
# Note that this script can be run together with Neptune logger
# In that case, user should specify the project and api-token in the train_model.py script
# If run together with neptune, subfolder will be created in the folder where model is stored

python3 src/train_model.py \
    # DATA: path to the expression data
    '/path/to/project/data/GTEx_GO_trimmed_expr.npy' \
    # PARAMETERS: path to json file containing model parameters (see example file)
    '/path/to/project/params/model_params.json' \
    # MASKS: path to masks indicating decoder connections
    '/path/to/project/ontologies/GO/GO_ensembl_trimmed_decoder_masks.pickle' \
    # MODELDIR: path where the best model will be stored
    '/path/to/project/models' \
    # LOG: 'true' if the run should be logged with Neptune, 'false' otherwise
    'true'




# This runs the script run_data_through_model.py
# Three output files would be generated from this call
# NAME_latent_space_embedding.npy - a 2d numpy array with samples in rows and latent layer terms in columns
# NAME_decoder_activities.npy - a 2d numpy array with samples in rows and terms in decoder in columns sorted by depth
# NAME_reconstruction.npy - a 2d numpy array with samples in rows and reconstructed genes in columns

python3 src/run_data_through_model.py \
    # DATA: path to the expression data
    '/path/to/project/data/GTEx_GO_trimmed_expr.npy' \
    # PARAMETERS: path to json file containing model parameters (see example file)
    '/path/to/project/params/model_params.json' \
    # MASKS: path to masks indicating decoder connections
    '/path/to/project/ontologies/GO/GO_ensembl_trimmed_decoder_masks.pickle' \
    # MODELPATH: path where the model is stored
    '/path/to/project/models/model' \
    # NAME: prefix to be used for saving files
    'GTEx'



# This runs the script log_images.py
# The script does the same as run_data_through_model.py but additionally logs images to the neptune run
# user should specify neptune project and api token in the script

python3 src/log_images.py \
    # DATA: path to the expression data
    '/path/to/project/data/GTEx_GO_trimmed_expr.npy' \
    # PARAMETERS: path to json file containing model parameters (see example file)
    '/path/to/project/params/model_params.json' \
    # MASKS: path to masks indicating decoder connections
    '/path/to/project/ontologies/GO/GO_ensembl_trimmed_decoder_masks.pickle' \
    # RUN: The neptune run
    'OVAE-113' \
    # MODELPATH: path where the model is stored
    '/path/to/project/models/model' \
    # ONTO_ANNOT: path to the annotation file of the ontology
    '/path/to/project/ontologies/GO/GO_ensembl_trimmed_annot.csv' \
    # SAMPLE_ANNOT: path to the annotation file for the samples
    '/path/to/project/data/GTEx_annot.csv' \
    # ANNOT_COLS: columns in sample to be used for plotting (up to two columns can be used, in that case should be separated by ,)
    'tissue' \
    # NAME: prefix to be used for saving files
    'GTEx'


