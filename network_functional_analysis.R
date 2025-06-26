# Author : Agatha Treveil
# Date : April 2020
#
# Functional overrepresentation analysis of causal networks: upstream signalling proteins (binding proteins-TFs inclusive), 
# and differentially expressed genes (DEGs) seperately
# Reactome and Gene ontology BP
#
# NB. GO analyses carried out using uniprot IDs, but for the Reactome analysis it is necessary to convert to ENTREZ IDs
#     All nodes of the contextualised specific PPI network (expressed omnipath) are used as the background for the upstream signalling proteins
#     All TGs of the contextualised specific TF-TG network (expressed dorothea) are used as the background for the DEGs
#     significantly overrepresented functions have q val <= 0.05
#     GO analyses use simplify to remove rudundant function (with default parameter 0.7)
#
# Input: whole network node file (output from combined_edge_node_tables.R)
#        background network file for ppis (contextualised specific PPI network output from filter_network_expressed_genes.R)
#        background network file for degs (contextualised specific TF-TG network output from filter_network_expressed_genes.R)
#        
# Output: Data table and dot plot for each overrepresentation analysis (commented out map plot due to package updates causing error when running it)


##### Set up #####

# Set timeout
options(timeout=900)


# Install requried packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("ReactomePA", quietly = TRUE))
  BiocManager::install("ReactomePA", update = FALSE, ask = FALSE)
if (!requireNamespace("clusterProfiler", quietly = TRUE)) 
  BiocManager::install("clusterProfiler")
if (!requireNamespace("tidyverse", quietly = TRUE)) 
  install.packages("tidyverse")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) 
  BiocManager::install("org.Hs.eg.db")

#packages
library(tidyverse)
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
outdir = 'func_anal_out'

##### Enrichment analysis #####

go_overrep <- function(net, name, id,folder){
  # Carries out GO overrepresentation analysis.

  # Carry out normal GO gene overrepresentation analysis
  go1 <- enrichGO(gene = as.character(net$Uniprot), OrgDb = "org.Hs.eg.db", keyType = id, ont = 'BP', qvalueCutoff = 0.05)
  
  if (nrow(go1) > 0) {
    # Remove redundancy of GO terms. Cutoff refers to similarity
    go2 <- clusterProfiler::simplify(go1, by = "qvalue", select_fun=min)
    
    # Get as dataframe
    go2_df <- as.data.frame(go2)
    
    #if (nrow(go2_df) >1){
    #  # Get enrichment map
    #  map_ora <- emapplot(go2)
    #  # Save map
    #  filen <- file.path(folder,paste0(name,"_map_GO_overrep.pdf"))
    #  pdf(filen)
    #  print(map_ora)
    #  dev.off()
    #}
    
    # Get dot plot
    dot_plot <- dotplot(go2, showCategory=10, orderBy="qvalue", font.size = 18)
    # Save dot plot
    filep <- file.path(folder,paste0(name,"_dot_GO_overrep.png"))
    png(filep, width = 225, height = 225, units='mm', res = 300)
    print(dot_plot)
    dev.off()
   
  } else{
    go2_df <- ""
  }
  
  return(go2_df)
}


reactome_overrep <- function(net, name, id, folder){
  # Carries out reactome overrepresentation analysis.
  
  # Convert to entrez genes
  net_e <- bitr(net$Uniprot, fromType=id, toType='ENTREZID', OrgDb="org.Hs.eg.db")
  #back_nodes_e <- bitr(back_nodes$node, fromType=id, toType='ENTREZID', OrgDb="org.Hs.eg.db")
  
  # Carry out normal reactome gene overrepresentation analysis (can only use entrez :( )
  re1 <- enrichPathway(gene = as.character(net_e$ENTREZ), organism = "human", qvalueCutoff = 0.05)
  
  if (nrow(re1) > 0) {
   
    # Get as dataframe
    re1_df <- as.data.frame(re1)
    
   #if (nrow(re1_df) >1){
   #  # Get enrichment map
   #  map_ora <- emapplot(re1)
   #  # Save map
   #  filen <- file.path(folder,paste0(name,"_map_reactome_overrep.pdf"))
   #  pdf(filen)
   #  print(map_ora)
   #  dev.off()
   #}
    
    # Get dot plot
    dot_plot <- dotplot(re1, showCategory=10, orderBy="qvalue", font.size = 10)
    # Save dot plot
    filep <- file.path(folder,paste0(name, "_dot_reactome_overrep.pdf"))
    pdf(filep, width=6.5)
    print(dot_plot)
    dev.off()
  } else{
    re1_df <- ""  
  }
  
  return(re1_df)
}


##### Run all #####

# Read nodes file
files<-list.files(path = "nodefiles", pattern=".csv")

for (i in files){
  nodes <- read.csv(file.path('nodefiles/',i), sep = ",", colClasses = c("NULL", "NULL", "NULL", "NULL", "NULL", "character"))
  part = strsplit(i, '_')
  part <- part[[1]][1]
  
  # Run GO overenrichment analysis
  go_res_ppi <- go_overrep(nodes, part, "UNIPROT", outdir)
  #if (go_res_ppi != ""){
  #  write.table(go_res_ppi, file = file.path(outdir,"go_overrep_results.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
  #}
  
  # Run reactome overenrichment analysis
  reactome_res_ppi <- reactome_overrep(nodes, part,"UNIPROT",outdir)
  #if (reactome_res_ppi != ""){
  #  write.table(reactome_res_ppi, file = file.path(outdir, "reactome_overrep_results.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
  #}
}
# reset message sink and close the file connection
sink(type="message")
close(zz)
