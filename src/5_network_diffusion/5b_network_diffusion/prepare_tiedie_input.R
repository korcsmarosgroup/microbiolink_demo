# Author: Agatha Treveil, Lejla Gul
# Date: June 2020 (modified in June, 2022)
#
# Script to prepare the input files for TieDIE
#
# Input: Contextualised PPI network output from 'filter_network_expressed_genes.R'
#        TF-DEGs network file output from 'get_regulator_deg_network.R'
#        Differentially expressed gene file output from 'diff_expression_deseq2.R'
#        Viral protein- human binding protein file (provided - based on Gordon et al.)
#
# Output: pathway.sif - sif format contextualised PPI (Omnipath) network
#         upstream.input - text file containing the upstream genes with a weight and direction (human binding partners)
#         downstream.intput - text file containing the downstream genes with a weight and direction (TFs)

##### Set up #####

# Install required packages
if (!requireNamespace("tidyverse", quietly = TRUE)) 
install.packages("tidyverse")

# Load packages
library(tidyverse)

outdir <- 'microbiolink_demo/results/TieDie/'


# Contextualised PPI network (omnipath) output from 'filter_network_expressed_genes.R'
ppis <- read.csv('microbiolink_demo/results/TieDie/1_process_a_priori_networks/omnipath_contextualised_network.txt', sep = "\t")

# TF-DEG interactions output from 'get_regulator_deg_network.R'
tfs <- read.csv('microbiolink_demo/results/TieDie/1_process_a_priori_networks/contextualised_regulator-deg_network.txt', sep = "\t")

# Filtered differential expression file output from 'diff_expression_deseq2.R'
degs <- read.csv('microbiolink_demo/results/TieDie/1_process_a_priori_networks/ec_ligands_epithelial_healthy.csv')

# Human binding proteins of bacterial proteins
hbps <- read.csv('microbiolink_demo/results/HMI/healthy_PPIs_all.csv', sep='\t', row.names = NULL)

# Create output dir if required 
path <- file.path(outdir, "2_network_diffusion", "input_files")
dir.create(path, showWarnings = FALSE, recursive=TRUE)


# Writinig out input files for 
hbps <- data.frame('human_protein' = hbps$row.names, 'bacterial_protein' = hbps$Bacterial.protein, 'sign' = '+')
write.table(hbps, file = file.path(path,"healthy_signed_PPIs.tsv"), sep = "\t",quote = FALSE, row.names = FALSE)

bac_annotation <- data.frame(unique(hbps$bacterial_protein),unique(hbps$bacterial_protein))
colnames(bac_annotation) <- c('uniprot', 'annotation')
write.table(bac_annotation, file = file.path(path,"bacterial_protein_annotation_healthy.csv"), sep = ",",quote = FALSE, row.names = FALSE)

##### Process pathways #####

# Convert to sif format
ppis2 <- ppis %>% mutate(direction = ifelse(consensus_stimulation == "1", "stimulates>", "inhibits>")) %>%
  dplyr::select(c(from,direction, to)) %>%
  unique()

# Save pathways
write.table(ppis2, file = file.path(path,"pathway.sif"), sep = "\t", col.names = F, row.names = F, quote = F)

##### Process upstream input #####

# Number of bacterial proteins bound to each human protein is the weight
# If sign of interaction given then use it, else assume all inhibitory ("-")
if ("sign" %in% colnames(hbps)){
  # Check values in sign column are "-" or "+"
  hbps_f <- hbps %>% filter((sign == "-") | (sign == "+"))
  if (nrow(hbps_f) != nrow(hbps)){
    print("WARNING: Some of the bacterial-human binding protein interactions were disgarded as the values in the 'sign' column were not '+' or '-'.")
  }
  if (nrow(hbps_f) == 0){
    stop("ERROR: bacterial-human binding protein interactions do not have the correct values in 'sign' column. They should be '+' or '-'.")
  }
  hbps2 <- hbps %>% dplyr::select(human_protein, sign) %>% dplyr::rename(direction=sign) %>% group_by(human_protein,direction) %>% summarise(n = n()) %>%
    dplyr::select(human_protein, n, direction)
} else {
  hbps2 <- hbps %>% dplyr::select(human_protein) %>% group_by(human_protein) %>% summarise(n = n()) %>% mutate(direction = "-")
}

# Save upstream data
write.table(hbps2, file = file.path(path,"upstream.input"), sep = "\t", col.names = F, row.names = F, quote = F)

##### Process downstream input #####

# (1/#targetgenes) * sum(lfc(targetgene)*signofint)

# Join the tf-deg network with the deg lfc values
tfs2 <- left_join(tfs, degs, by =c("to"="Gene")) %>% dplyr::select(c(from, to, consensus_stimulation, Expression)) %>%
  mutate(lfc_sign = ifelse(consensus_stimulation == "1", Expression, -Expression))

# Get the sum of all lfc*sign values - and the number of target genes for each tf
tfs3 <- tfs2 %>% dplyr::select(from, lfc_sign) %>% group_by(from) %>% summarise(sumof = sum(lfc_sign), n = n())

# Divide sumof by n and determine sign (based on sign of the value)
tfs4 <- tfs3 %>% mutate(final_val = sumof/n) %>% mutate(sign = ifelse((final_val >= 0), "+", "-")) %>%
  dplyr::select(from, final_val, sign)

# Save downstream data
write.table(tfs4, file = file.path(path,"downstream.input"), sep = "\t", col.names = F, row.names = F, quote = F)
