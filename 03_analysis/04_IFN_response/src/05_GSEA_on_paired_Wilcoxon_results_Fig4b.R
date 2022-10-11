
# import libraries and set dirs

library(tidyverse)
library(clusterProfiler)


project_dir <- getwd()


# load results from paired Wilcoxon testing and order by pvalue
wilcox <- read.csv(paste0(project_dir, '/03_analysis/04_IFN_response/results/Genes_upregulating_IFN_activity_Wilcoxon_results.csv'))

# load genes upregulated in IFN-ß treated CD4T cells
diff_genes <- read.csv(paste0(project_dir, '/03_analysis/04_IFN_response/results/CD4T_IFN-ß_stim_up.txt'), header=FALSE)

# prepare input for clusterProfiler
term2gene <- data.frame(term = 'CD4T_stim_up', gene = diff_genes)
genelist <- setNames(-log10(wilcox$w_pval), wilcox$gene)

# run GSEA and save the results
res <- clusterProfiler::GSEA(geneList=genelist, TERM2GENE = term2gene)
saveRDS(res, paste0(project_dir, '/03_analysis/04_IFN_response/results/reference_genes_predicted_genes_GSEA.RDS'))

# make GSEA plot and save
p <- enrichplot::gseaplot2(res, geneSetID = 1, pvalue_table = TRUE, ES_geom = "line")
ggsave(paste0(project_dir, '/03_analysis/04_IFN_response/figures/05_IFN_signaling_Wilcoxon_GSEA_Fig3b.pdf'), p)

# leading edge cutoff: 717


