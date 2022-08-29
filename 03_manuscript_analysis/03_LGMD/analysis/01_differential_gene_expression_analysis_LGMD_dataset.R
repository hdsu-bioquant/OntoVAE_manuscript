
library(tidyverse)
library(DESeq2)

project_dir <- '/net/data.isilon/ag-cherrmann/ddoncevic/projects/OntoVAE_manuscript'


# load datas
expr <- read.table(paste0(project_dir, '/datsets/LGMD/GSE202745_exp3105-x2-switched-RawCountsTable.txt'), sep="\t", header=TRUE)
expr_sub <- expr[,c(2:86)]


# remove lowly expressed features
rs <- rowSums(expr_sub[,-1])
expr_sub <- expr_sub[rs > 50,]

# remove duplicate gene names
dups <- expr_sub$Gene.Name[duplicated(expr_sub$Gene.Name)]
expr_sub <- expr_sub[!expr_sub$Gene.Name %in% dups,]
rownames(expr_sub) <- expr_sub$Gene.Name
expr_sub <- expr_sub[,-1]

# create annotation table
annot <- data.frame(ID=colnames(expr_sub)) %>% separate(ID, c('Status', 'Number', 'Tissue'))
rownames(annot) <- colnames(expr_sub)

# remove the PCA outliers
outliers <- c("C_9_SM", "P_3_SM", "P_1_VLM", "P_5_SM", "P_5_VLM", "P_10_SM", "P_12_SM", "P_12_VLM", "P_15_SM")

# generate DESeq dataset
dds <- DESeqDataSetFromMatrix(countData=expr_sub[,!colnames(expr_sub) %in% outliers],
                             colData=annot[!rownames(annot) %in% outliers,],
                             design= ~Status)

# get differentially expressed genes
dds <- DESeq(dds)
res <- lfcShrink(dds, contrast=c("Status", "P", "C"), alpha=0.05, type="ashr")

# write results to file
write.csv(res, paste0(project_dir, '/03_manuscript_analysis/03_LGMD/results/LGMD_DGA_results.csv'), quote=FALSE)


