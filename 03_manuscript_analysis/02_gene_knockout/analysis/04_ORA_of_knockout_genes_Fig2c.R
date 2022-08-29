
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(ggplot2)


working_dir <- "/net/data.isilon/ag-cherrmann/ddoncevic/projects/OntoVAE_manuscript"


###----------------> define experiment genes

genes <- c('DMD', 'SFSWAP', 'COX5A')


###----------------> import Wilcoxon results

res <- read.csv(paste0(working_dir, "/03_manuscript_analysis/02_gene_knockout/results/GTEx_muscle_OVAE-113_gene_knockout_reconstruction_Wilcoxon_results.csv"), sep=";")


###---------------->  function to convert between different gene ids

id_convert <- function(id_list, from, to){
    if(length(id_list) > 0){
        gene <- select(org.Hs.eg.db, id_list, c(to), from)
        return(gene[to][!is.na(gene[to])])
    } else {
        return(NULL)
    }
}


###----------------> function to get genelist

get_genelist <- function(res, gene){
    res <- res[res$ko_gene == gene,]
    res <- res[order(res$qval),]
    genelist <- res$gene[1:100]
    return(id_convert(as.character(genelist), 
                        "ENSEMBL", 
                        "ENTREZID"))
}




###----------------> create genelists for GO enrichment 

gene_lists <- lapply(genes, function(g) get_genelist(res, g))
names(gene_lists) <- genes


###----------------> helper function for plotting

overrep_bar_plot <- function(res) {

    res <- data.frame(res)

    col <- 'darkgrey'

    if(nrow(res) > 10){
        res <- res[1:10,]
    }
    res$Description <- factor(res$Description, levels = rev(res$Description))

    p <- ggplot(res, aes(x=-log10(pvalue), y=Description)) + 
         geom_bar(stat='identity', fill=col, width=0.5) +
         labs(y="") +
         theme_classic() +
         geom_vline(xintercept=-log10(0.05), color='darkred', linetype='dashed') +
         theme(axis.text=element_text(size=12))
        
    return(p)
}


###----------------> function to make the plots


ora_plots <- function(genelist, gene, name){
    res_go <- enrichGO(gene = genelist,
                    OrgDb = 'org.Hs.eg.db', 
                    ont = "BP", 
                    pvalueCutoff = 0.05)
    #write.csv2(data.frame(res_go), paste0(working_dir, "/03_manuscript_analysis/02_gene_knockout/results/GO_ORA_reconstruction_after_", gene, "_knockout.csv"), row.names=FALSE, quote = FALSE)
    go_bar <- overrep_bar_plot(res_go)
    ggsave(paste0(working_dir, "/03_manuscript_analysis/02_gene_knockout/figures/GO_ORA_reconstruction_after_", gene, "_knockout.pdf"), go_bar, width=15, height=8)
}



# run the functions to create the plots
for(g in genes) {
     print(g)
     ora_plots(gene_lists[[g]], g, "")
 }



