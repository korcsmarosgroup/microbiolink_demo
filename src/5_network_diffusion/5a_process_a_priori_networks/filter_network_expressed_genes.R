# Author: Agatha Treveil, Lejla Gul
# Date: April 2020 (modified in June, 2022)
#
# Script to filter networks for expressed genes - runs for regulatory interactions and for ppis
#
# Input: Network file, space delimited, gene names in columns "source_genesymbol" "target_genesymbol" or "to" "from"
#        Table of genes which are expressed, tab delimited, gene names in column "Gene" (ouput from DESeq2)
#        ID type of the differentially expressed genes - uniprot or gene symbols
#
# Output: Networks in same format as input network (but tab seperated), filtered to include only interactions where source and target node are in the expressed list.

# Installing packages
if (!requireNamespace("tidyverse", quietly = TRUE)) 
  install.packages("tidyverse")

# Load packages
library(tidyverse)


# Output directory
outdir <- 'microbiolink_demo/results/TieDie/'

# ID type of the  expressed genes - uniprot or gene symbols
id_type <- 'symbol' # symbol or uniprot - the ids in the expression data

# Create output dir if required
path <- file.path(outdir, "1_process_a_priori_networks")
dir.create(path, showWarnings = FALSE, recursive=TRUE)

if (id_type == "symbol"){
  source_col = "source_genesymbol"
  target_col = "target_genesymbol"
} else if (id_type == "uniprot") {
  source_col = "to"
  target_col = "from"
}

# Gene expression file - tab delimited
expressed <- read.csv('microbiolink_demo/results/human_proteins/healthy_expressed_genes.csv', sep = ",")
dorothea <- 'microbiolink_demo/results/TieDie/1_process_a_priori_networks/unprocessed_networks/dorothea_abc_signed_directed.txt'
omnipath <- 'microbiolink_demo/results/TieDie/1_process_a_priori_networks/unprocessed_networks/omnipath_signed_directed.txt'

files <- c(dorothea, omnipath)
  
for (i in files){
  
  # Network file - space delimited 
  network <- read.csv(file.path(i), sep = " ")
  
  # Filter source and target nodes
  network_f <- network %>% filter((get(source_col) %in% expressed$Gene) & (get(target_col) %in% expressed$Gene))
  
  # Get network name for out filename - depends on the file path
  file <- strsplit(i, "/")[[1]][13]
  name <- strsplit(file, "_")[[1]][1]
  
  # Save output
  write.table(network_f, file = file.path(path, paste0(name, "_contextualised_network.txt")), sep = "\t", quote = F, row.names = F)
}

