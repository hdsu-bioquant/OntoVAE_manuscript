library(tidyverse)

project_dir <- getwd()


expr <- read.table(paste0(project_dir, '/datasets/LGMD/GSE202745_exp3105-x2-switched-RawCountsTable.txt'), sep="\t", header=TRUE)

tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

expr_tpm <- tpm(expr[,c(3:86)], expr$Length)

log_tpm <- log2(expr_tpm + 1)
rownames(log_tpm) <- expr$Gene.ID


annot <- data.frame(ID=colnames(log_tpm)) %>% separate(ID, c('Status', 'Number', 'Tissue'))


pca <- prcomp(log_tpm, center=TRUE, scale = TRUE)
pca_res <- pca$rotation[,c(1:2)]

plot_in <- data.frame(cbind(pca_res, annot))

p <- ggplot(plot_in, aes(x=PC1, y=PC2, color=colSums(log_tpm), shape=Status)) + 
     geom_point() 

ggsave(paste0(project_dir, '/03_anayĺysis/03_LGMD/figures/01_PCA_on_LGMD_dataset.png'))


outliers <- rownames(plot_in[plot_in$PC1 > -0.1 & plot_in$PC2 > 0,])

#outliers: "C_9_SM", "P_3_SM", "P_1_VLM", "P_5_SM", "P_5_VLM", "P_10_SM", "P_12_SM", "P_12_VLM", "P_15_SM"