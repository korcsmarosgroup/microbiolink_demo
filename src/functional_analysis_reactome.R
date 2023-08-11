# Author: Agatha Treveil, Lejla Gul
# Updated: May 2023
# Functional analysis using Reactome DB
# Input: output from MicrobioLink (list of bacteria targeted proteins) - target
#        output from Location analysis from OmniPath (list of membrane-based cell specific proteins) - background
#        name of the output dotplot PDF
#        name of the output barplot RDS
#
# Output: doptlpot PDF
#         barplot RDS

##### Set up #####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ReactomePA")

BiocManager::install("clusterProfiler")

BiocManager::install("AnnotationDbi")

#Human data
BiocManager::install("org.Hs.eg.db") 

install.packages("ggplot2")

#packages
library(dplyr)
library(ReactomePA)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)


##### Overrepresentation analysis #####

# Empty df's for data
all_reactome <- data.frame()

overrep <- function(genes, universe_file, dotplot, barplot){
  # Remove duplicates
  genes2 <- genes %>% unique() %>% na.omit()
  
  # Convert to gene IDs if it is necessary
  genes2_e <- bitr(genes2$Human.UniProt, fromType='UNIPROT', toType='ENTREZID', OrgDb="org.Hs.eg.db")
  
  #### REACTOME #### 
  
  #Set up the background gene set - all the expressed membrane proteins
  universe_table <- read.csv(universe_file)
  universe_table
  universe2 <- universe_table %>% unique() %>% na.omit()
  universe_e <- bitr(universe2$Human.Gene.Symbol, fromType='SYMBOL', toType='ENTREZID', OrgDb="org.Hs.eg.db")
  
  
  #Reactome pathway enrichment
  genes_reactome <- enrichPathway(gene=genes2_e$ENTREZID,organism='human',qvalueCutoff = 0.5, pvalueCutoff = 0.5, readable=T, universe = universe_e$ENTREZID)
  
  #Get reactome results
  genes_reactome2 <- as.data.frame(genes_reactome)
  #genes_reactome2 <- mutate(genes_reactome2, group = name)
  
  #Append reactome results to previous results
  all_reactome <- rbind(all_reactome, genes_reactome2)
  
  # Get dot plot
  dot_plot_r <- dotplot(genes_reactome, showCategory=20, font.size = 8)
  
  # Save dotplot
  pdf('Reactome_EEC_H_HMI_targets.pdf')
  print(dot_plot_r)
  dev.off()
  
  barplot_r <- barplot('barplot_Reactome_EEC_H_HMI_targets.rds', x="qvalue")
  
  # Optional changing the visual style of the barplot
  #barplot(genes_reactome_allgene) + scale_fill_viridis() + theme_classic()
  
  saveRDS(barplot_r, file=barplot)

  return(all_reactome)
}

all_results_all_lists <- data.frame()

# Run the function
genes<- read.csv("list of genes", header=TRUE)
all_reactome <- overrep(genes, 'background', "output_dotplot.pdf", 'output_barplot.rds')
universe_table <- read.csv("background_geneset.csv")
write.table(all_reactome, "Reactome_output.txt", sep = "\t", row.names = FALSE, quote = FALSE)

