library(tidyverse)
library(ComplexHeatmap)
library(viridis)
library(RColorBrewer)

project_dir <- paste0(getwd(), '/03_analysis/01_pathway_activities')

# load data
aucs <- read.csv(paste0(project_dir, "/results/pathway_activities_tissue_classification_aucs.csv"), sep=";")
onto_annot <- read.csv(paste0(project_dir, "/data/GO_annot.csv"), sep=";")

# get top 10 terms per tissue
top_terms <- aucs %>% 
    pivot_longer(-Name) %>%
    group_by(name) %>% 
    slice_max(order_by = value, n = 10) %>%
    filter(value > 0.5)
top_terms <- unique(top_terms$Name)

# filter auc matrix for top terms
top_aucs <- aucs[aucs$Name %in% top_terms,]
rownames(top_aucs) <- top_aucs$Name
top_aucs <- top_aucs[,-1]

# filter onto annot
annot_filt <- onto_annot[onto_annot$Name %in% top_terms,]

# define heatmap colors
heat_col = colorRamp2(quantile(heat_data, seq(0,1,by=0.12)), rev(brewer.pal(9, "RdYlBu")))  
my_cols <- list(Depth = setNames(colorRampPalette(brewer.pal(9, "Greys"))(length(unique(annot_filt$depth))), unique(annot_filt$depth)))
  
# define annotations
col_annot <- HeatmapAnnotation(Tissue = colnames(top_aucs)), col=my_cols)
row_annot <- rowAnnotation(Depth = annot_filt$depth, col=my_cols)

# plot the heatmap
heatmap <- Heatmap(top_aucs,
                   name="AUC",
                   col=magma(n=9),
                   column_title=NULL,
                   cluster_columns=FALSE, 
                   #show_row_dend=FALSE,
                   show_row_names=FALSE,
                   show_column_names=FALSE,
                   border=TRUE,
                   top_annotation=col_annot,
                   right_annotation=row_annot)

pdf(paste0(project_dir, "/figures/04_Top_AUC_heatmap.pdf"))
heatmap
dev.off()