while read p; do
    echo "$p"
    python3 src/10_perform_systematic_knockout_in_models_with_removed_link.py \
        "$p"
done <results/LGMD_genes_leading_edge.txt