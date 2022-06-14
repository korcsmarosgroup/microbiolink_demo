# Author : Agatha Treveil, Lejla Gul
# Date : April 2020 (modified in June 2022)
#
# Script to combine subnetworks into 1 network table and one node table
#
# Input: Receptor-tf network output from tiedie
#        Bacteria - human binding protein network (columns 'bacterial_protein', 'human_protein')
#        Contextualised regulatory network (TFs - DEGs) (incl. columns 'DEG',"	TF", "consensus_stimulation") output from 'get_regulator_deg_network.R'
#        Heats values output from TieDie
#        Bacterial proteins gene symbol to uniprot conversion table
#        Differential expression table (unfiltered)
#
# Output: Network file where each line represents an interaction
#         Node table where each lines represents one node in the network

##### Setup #####

# Set timeout
options(timeout=120)

# Capture  messages and errors to a file.
# Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("tidyverse", quietly = TRUE)) 
  install.packages("tidyverse")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")

library(tidyverse)
library(org.Hs.eg.db)

# Define parameters

outdir <- 'microbiolink_demo/results/TieDie/'

# Receptor-tf network from tiedie
rec_tf <- read.csv('microbiolink_demo/results/TieDie/tiedie.cn.sif', header=F, col.names=c("Source.node","Relationship","Target.node"), sep = "\t")

# Node heats from tiedie
heats <- read.csv('microbiolink_demo/results/TieDie/heats.NA', sep="=")

# Bacterial protein-receptor network
hbps <- read.csv('microbiolink_demo/results/TieDie/2_network_diffusion/input_files/healthy_signed_PPIs.tsv', sep = "\t")

# Bacterial proteins Uniprot IDs
bacteria <- read.csv('microbiolink_demo/results/TieDie/2_network_diffusion/input_files/bacterial_protein_annotation_healthy.csv', sep = ",")

# tf-deg network
tf_deg <- read.csv('microbiolink_demo/results/TieDie/1_process_a_priori_networks/contextualised_regulator-deg_network.txt', sep = "\t")

# expression table (unfiltered)
exp <- read.csv('microbiolink_demo/results/human_proteins/healthy_expressed_genes.csv')

# ID type of the expression table (uniprot or symbol)
id_type <- 'symbol' # "symbol" or "uniprot"

# Create output dir if required
path <- file.path(outdir, "3_create_network")
dir.create(path, showWarnings = FALSE, recursive=TRUE)

##### Process network edges #####

rec_tf2 <- rec_tf %>% mutate(layer= "bindingprot-tf") 

if("sign" %in% colnames(hbps)){
  hbps2 <- hbps %>% mutate(Relationship = ifelse(sign == "-","inhibits>", ifelse(sign=="+", "stimulates>", "unknown")))  %>%
    dplyr::select(-c(sign)) %>%
    mutate(layer = "bacteria-bindingprot") %>%
    dplyr::rename(Source.node = bacterial_protein, Target.node = human_protein) %>%
    filter(Target.node %in% rec_tf2$Source.node)
} else {
  hbps2 <- hbps %>% mutate(Relationship = "unknown", layer = "bacteria-bindingprot") %>%
    dplyr::rename(Source.node = bacterial_protein, Target.node = human_protein) %>%
    filter(Target.node %in% rec_tf2$Source.node)
}

tf_deg2 <- tf_deg %>% mutate(layer = "tf-deg") %>%
  dplyr::select(Target.node = to, Source.node = from, Relationship = consensus_stimulation, layer) %>%
  filter(Source.node %in% rec_tf2$Target.node) %>%
  mutate(Relationship = str_replace(Relationship, "1", "stimulates>")) %>%
  mutate(Relationship = str_replace(Relationship, "0", "inhibits>"))

# Join together
whole_net <- rbind(hbps2, rec_tf2, tf_deg2) %>% unique()

# Save
write.table(whole_net, file = file.path(path, "final_network.txt"), sep = "\t", quote = F, row.names = F)

##### Process node table #####

# Get layers of all nodes
nodes1 <- whole_net %>% filter(layer == "bacteria-bindingprot") %>% mutate(bacteria_layer = "bacteria") %>% dplyr::select(node = Source.node, bacteria_layer) %>% unique()
nodes2 <- whole_net %>% filter(layer == "bacteria-bindingprot") %>% mutate(bindingprot_layer = "bindingprot") %>% dplyr::select(node = Target.node, bindingprot_layer) %>% unique()
nodes3 <- whole_net %>% filter(layer == "bindingprot-tf") %>% mutate(ppi1_layer = "bindingprot and/or protein") %>% dplyr::select(node = Source.node, ppi1_layer) %>% unique()
nodes4 <- whole_net %>% filter(layer == "bindingprot-tf") %>% mutate(ppi2_layer = "protein and/or tf") %>% dplyr::select(node = Target.node, ppi2_layer) %>% unique()
nodes5 <- whole_net %>% filter(layer == "tf-deg") %>% mutate(tf_layer = "tf") %>% dplyr::select(node = Source.node, tf_layer) %>% unique()
nodes6 <- whole_net %>% filter(layer == "tf-deg") %>% mutate(deg_layer = "deg") %>% dplyr::select(node = Target.node, deg_layer) %>% unique()

# Join node layers together
all_nodes <- full_join(nodes1, nodes2) %>% full_join(nodes3) %>% full_join(nodes4) %>% full_join(nodes5) %>% full_join(nodes6)

# Combine the ppi layer col
all_nodes <- all_nodes %>% mutate(ppi_layer = ifelse((ppi1_layer == "bindingprot and/or protein" & ppi2_layer == "protein and/or tf"), "protein", "NA")) %>%
  dplyr::select(-c(ppi1_layer, ppi2_layer))

# Combine into 1 column
all_nodes <- all_nodes %>% unite(all_nodes, bacteria_layer, bindingprot_layer, ppi_layer, tf_layer, deg_layer, sep = ";", remove=FALSE, na.rm = TRUE)
all_nodes[is.na(all_nodes)] <- "NA"
rm(nodes1,nodes2, nodes3, nodes4, nodes5, nodes6)

# Get human node id conversion from Orgdb
id_mapping <- select(org.Hs.eg.db, keys=all_nodes$node, columns=c("SYMBOL","ENTREZID"), keytype="UNIPROT")
# Get only the first mapped id - as need 1:1 mapping later on
id_mapping2 <- id_mapping[!(duplicated(id_mapping$UNIPROT)),]
# Add id conversions to node table
all_nodes <- left_join(all_nodes, id_mapping2, by = c("node"="UNIPROT"))

# Add bacteria node id conversions
bacteria_sym <- bacteria %>%dplyr::select(c(node = annotation, uniprot))
all_nodes <- left_join(all_nodes, bacteria_sym, by = c("node")) %>% mutate(SYMBOL = replace_na(SYMBOL,"")) %>% mutate(uniprot = replace_na(as.character(uniprot),""))
all_nodes <-all_nodes %>% mutate(SYMBOL = paste0(SYMBOL, uniprot)) %>% dplyr::select(-c(uniprot))

# Add node heats from tiedie
heats$node <- row.names(heats)
heats$node <- gsub('\\s+', '', heats$node)
all_nodes <- left_join(all_nodes,heats)

# Add expression values
if ("Gene" %in% colnames(exp)){
  exp <- exp %>% dplyr::rename(X = Gene)
}

if (id_type == "symbol"){
  exp2 <- exp %>% dplyr::select(SYMBOL = X, Expression)
} else if(id_type == "uniprot"){
  exp2 <- exp %>% dplyr::select(node = X, Expression)
}
all_nodes <- left_join(all_nodes, exp2)

# Save node table
write.table(all_nodes, file = file.path(path, "node_table.txt"), sep = "\t", quote = F, row.names = F)


