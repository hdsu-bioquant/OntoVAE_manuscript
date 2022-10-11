while read p; do
    echo "$p"
    python3 src/09_train_GTEx_HPO_models_with_removed_link.py \
        'HP:0006785' \ # Limb-girdle muscular dystrophy
        "$p"
done <results/LGMD_genes_leading_edge.txt