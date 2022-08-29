

library(tidyverse)

working_dir <- '/net/data.isilon/ag-cherrmann/ddoncevic/projects/OntoVAE_manuscript'


#load muscle experiment genes

genes <- c('DMD', 'SFSWAP', 'COX5A')


# load Wilcoxon results

res <- read.csv(paste0(working_dir, '/03_manuscript_analysis/02_gene_knockout/results/GTEx_muscle_OVAE-113_gene_knockout_GO_decoder_Wilcoxon_results.csv'), sep=";")



# function to generate bar plot of top 10 terms

overrep_bar_plot <- function(gene) {

    data <- res[res$gene == gene,] 
    data <- data[order(data$qval),]

    data_sub <- data[1:10,]

    col <- 'darkgrey'

    data_sub$term <- factor(data_sub$term, levels = rev(data_sub$term))

    p <- ggplot(data_sub, aes(x=log10p, y=term)) + 
         geom_bar(stat='identity', fill=col, width=0.5) +
         labs(y="") +
         theme_classic() +
         geom_vline(xintercept=-log10(0.05), color='darkred', linetype='dashed') +
         theme(axis.text=element_text(size=12))
        
    return(p)
}


# loop through genes and generate ggplots

for(gene in genes){
    p <- overrep_bar_plot(gene)
    ggsave(paste0(working_dir ,'/03_manuscript_analysis/02_gene_knockout/figures/GO_terms_downregulated_after_', gene, '_knockout.pdf'), p,  width=15, height=8)
}