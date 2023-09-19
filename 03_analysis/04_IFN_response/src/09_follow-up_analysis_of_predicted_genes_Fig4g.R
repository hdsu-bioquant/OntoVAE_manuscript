
library(tidyverse)
library(rjson)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(RColorBrewer)
library(wesanderson)

project_dir <- getwd()
msig_dir <- '/net/data.isilon/ag-cherrmann/ddoncevic/data/MSigDB' 

# function for id conversion

id_convert <- function(id_list, from, to){
    if(length(id_list) > 0){
        gene <- select(org.Hs.eg.db, id_list, c(to), from)
        return(gene[to][!is.na(gene[to])])
    } else {
        return(NULL)
    }
}


# import the gene lists

res <- read.csv(paste0(project_dir, '/03_analysis/04_IFN_response/results/Genes_upregulating_IFN_activity_Wilcoxon_results.csv'))
predict_genes <- res$gene[1:717]

diff_genes <- read.table(paste0(project_dir, '/03_analysis/04_IFN_response/results/CD4T_IFN-ß_stim_up.txt'))
pred_genes <- predict_genes[!predict_genes %in% diff_genes$V1]


# load MSigDB datasets

c7 <- read.gmt(paste0(msig_dir, "/c7.all.v7.5.1.symbols.gmt"))


# grep for Interferon terms in C7

inf_c7 <- c7[grep("IFN", c7$ont),]

# perform overrepresentation analysis

res_c7 <- enricher(pred_genes, TERM2GENE=inf_c7)


# make the dotplot
pal <- rev(wes_palette("Zissou1", 100, type = "continuous"))
p <- dotplot(res_c7) + scale_colour_gradientn(colours = pal)

# save the dotplot
ggsave(paste0(project_dir, '/03_analysis/04_IFN_response/figures/08-01_Enrichment_of_IFN_related_genesets.pdf'), p, width=11, height=5)





