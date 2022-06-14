# Author: Agatha Treveil, Lejla Gul
# Date: May 2020 (modified in June 2022)
#
# Script to get regulator - differentially expressed genes from the contextualised regulatory network
#
# Input: Contextualised regulatory network file (all nodes expressed, tab delimited, output from filter_network_expressed_genes.R)
#           uniprot ids in columns "to" and "from", gene symbols in columns "source_genesymbol" and "target_genesymbol"
#        Table of differentially expressed genes (csv, output from deseq2)
#        ID type of the differentially expressed genes - uniprot or gene symbols
#
# Output: Network in same format as input network (but tab seperated), filtered to include only interactions where target nodes are differentially expressed.

##### Set up #####

# Install required packages
if (!requireNamespace("tidyverse", quietly = TRUE)) 
  install.packages("tidyverse")

# Load required packages
library(tidyverse)

# Output directory
outdir <- 'microbiolink_demo/results/TieDie/'

# Create output dir if required
path <- file.path(outdir, "1_process_a_priori_networks")
dir.create(path, showWarnings = FALSE, recursive=TRUE)

# Contextualised regulatory network
reg_net <- read.csv('microbiolink_demo/results/TieDie/1_process_a_priori_networks/omnipath_contextualised_network.txt', sep = "\t")
    
# Expressed genes (prefiltered)
endpoint_genes <- read.csv('microbiolink_demo/results/TieDie/1_process_a_priori_networks/ec_ligands_epithelial_healthy.csv')

# ID type of the differentially expressed genes - uniprot or gene symbols
id_type <- 'uniprot' # either "uniprot" or "symbol"
  
##### Preprocess #####

# Get column names of network to match to the differentially expressed genes (based on id type)
if(id_type == "symbol"){
  source_col <- "source_genesymbol"
  target_col <- "target_genesymbol"
} else if(id_type == "uniprot"){
  source_col <- "from"
  target_col <- "to"
} else {
  stop("The differential expression data id type is not correctly specified. Should be \"uniprot\" or \"symbol\"")
}

##### Filter network #####

# Filter netowrk so all target/endpoint genes are expressed
reg_net_f <- reg_net %>% filter(get(target_col) %in% endpoint_genes$Gene) %>% unique()

# Get list of regulators with number of targeted DEGs
regs <- reg_net_f %>% dplyr::select(deg_regs = !!source_col) %>% add_count(deg_regs, name = "num_degs") %>% unique()

##### Save #####

write.table(reg_net_f, file = file.path(path,"contextualised_regulator-deg_network.txt"), sep = "\t", quote = F, row.names = F)
write.table(regs, file = file.path(path, "contextualised_regulators_of_degs.txt"), sep = "\t", quote = F, row.names = F)
