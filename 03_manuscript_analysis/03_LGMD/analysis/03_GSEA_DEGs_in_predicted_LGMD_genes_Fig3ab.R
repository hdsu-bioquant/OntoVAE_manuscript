
library(tidyverse)
library(clusterProfiler)
library(rjson)

working_dir <- '/net/data.isilon/ag-cherrmann/ddoncevic/projects/OntoVAE_manuscript'


# DGA results

res <- read.csv(paste0(working_dir, '/03_manuscript_analysis/03_LGMD/results/LGMD_DGA_results.csv'))


# knockout results


ko_res_up <- read.csv(paste0(working_dir, '/03_manuscript_analysis/03_LGMD/results/GTEx_muscle_OVAE-120_gene_knockout_LGMD_Wilcoxon_results_up.csv'))
ko_res_dn <- read.csv(paste0(working_dir, '/03_manuscript_analysis/03_LGMD/results/GTEx_muscle_OVAE-120_gene_knockout_LGMD_Wilcoxon_results_dn.csv'))


# filter the DGA results

res_filt <- res[res$X %in% ko_res_up$gene,]

res_sig <- res_filt[!is.na(res_filt$padj),]
res_sig <- res_sig[res_sig$padj < 0.1,]

# split in LGMD_up and LGMD_dn

lgmd_up = res_sig$X[res_sig$log2FoldChange > 0]
lgmd_dn = res_sig$X[res_sig$log2FoldChange < 0]

write.table(lgmd_up, paste0(working_dir, '/03_manuscript_analysis/03_LGMD/results/LGMD_up.txt'), row.names=FALSE, col.names=FALSE,quote=FALSE)
write.table(lgmd_dn, paste0(working_dir, '/03_manuscript_analysis/03_LGMD/results/LGMD_dn.txt'), row.names=FALSE, col.names=FALSE,quote=FALSE)

# GSEA

term2gene <- rbind(data.frame(term = 'LGMD_up', gene = lgmd_up),
                   data.frame(term='LGMD_dn', gene = lgmd_dn))

genelist_up <- setNames(-log10(ko_res_up$w_qval), ko_res_up$gene)
genelist_dn <- setNames(-log10(ko_res_dn$w_qval), ko_res_dn$gene)


gsea_res_up <- clusterProfiler::GSEA(geneList=genelist_up, TERM2GENE = term2gene)
gsea_res_dn <- clusterProfiler::GSEA(geneList=genelist_dn, TERM2GENE = term2gene)

p <- enrichplot::gseaplot2(gsea_res_up, geneSetID = 1, pvalue_table = TRUE, ES_geom = "line")
ggsave(paste0(working_dir, '/03_manuscript_analysis/03_LGMD/figures/03-01_GSEA_HPO_LGMD_up.pdf'), p)

p <- enrichplot::gseaplot2(gsea_res_dn, geneSetID = 1, pvalue_table = TRUE, ES_geom = "line")
ggsave(paste0(working_dir, '/03_manuscript_analysis/03_LGMD/figures/03-01_GSEA_HPO_LGMD_dn.pdf'), p)



# write the lgmd genes to file

up_lgmd_genes <- Reduce(intersect, list(names(genelist_up[1:424]), lgmd_genes))
dn_lgmd_genes <- Reduce(intersect, list(names(genelist_dn[1:564]), lgmd_genes))


write.csv(data.frame(c(up_lgmd_genes, dn_lgmd_genes)), paste0(working_dir, '/03_manuscript_analysis/03_LGMD/results/LGMD_genes_in_leading_edge.csv'), row.names=FALSE, col.names=FALSE, quote=FALSE)

