# ===============================================================================
# Pathway enrichment analysis and functional module identification
# Xueqing Jia, 2025
# Description:
#  Identify and visualize biological functional modules from significant proteins
#  through GO, KEGG, and Reactome enrichment, similarity-based clustering,
#  and network-based community detection.
# ================================================================================

### Load R packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(DOSE)
library(enrichplot)
library(cowplot)
library(simplifyEnrichment)
library(reshape2)
library(igraph)
library(dplyr)
library(ReactomePA)
library(tibble)
library(tidyr)
library(purrr)

### Load significant proteins (Model 2)
sig_marker<-read.csv(file = "lm_proteins&FI_Model2.csv",row.names = 1)
module_gene<-subset(sig_marker,sig_marker$P_bon<=0.05)#FI_sum 1339

### Convert gene SYMBOL to ENTREZID
ids <- bitr(module_gene,fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db")
genes = ids[,2] # list of ENTREZ IDs

#################### GO enrichment ###########################
enrich_go <- enrichGO(gene=genes,'org.Hs.eg.db',
                      ont="ALL",
                      pAdjustMethod = "BH",
                      minGSSize = 10,
                      pvalueCutoff=0.05,
                      qvalueCutoff  = 0.05,
                      readable=TRUE)
enrich_go_ALL<-data.frame(enrich_go)

### Calculate similarity among GO terms
sim_BP <- GO_similarity(filter(enrich_go, ONTOLOGY == 'BP')$ID, ont = 'BP', db = 'org.Hs.eg.db', measure = "Wang")
sim_CC <- GO_similarity(filter(enrich_go, ONTOLOGY == 'CC')$ID, ont = 'CC', db = 'org.Hs.eg.db', measure = "Wang")
sim_MF <- GO_similarity(filter(enrich_go, ONTOLOGY == 'MF')$ID, ont = 'MF', db = 'org.Hs.eg.db', measure = "Wang")

## Convert similarity matrices to long format
sim_BP_corr.table <- melt(replace(sim_BP, lower.tri(sim_BP, TRUE), NA), na.rm = TRUE)
sim_CC_corr.table <- melt(replace(sim_CC, lower.tri(sim_CC, TRUE), NA), na.rm = TRUE)
sim_MF_corr.table <- melt(replace(sim_MF, lower.tri(sim_MF, TRUE), NA), na.rm = TRUE)

sim_ALL_corr.table <- rbind(sim_BP_corr.table,sim_CC_corr.table,sim_MF_corr.table)

### Select similarity score greater than 0.5 to construct GO term similarity network
data <- sim_ALL_corr.table %>%
  filter(value >= 0.5, Var1 != Var2)
head(data)

## Compute degree (number of connections) for each node
vertices<-c(as.character(data$Var1), as.character(data$Var2)) %>% as_tibble()
vertices<-data.frame(table(vertices$value))
colnames(vertices) <- c("node","n")

## Add GO ontology type information
type <- enrich_go_ALL[,c("ID","ONTOLOGY")]
colnames(type)<-c("node","type")

vertices %>%
  left_join(type, by="node")-> vertices

## Build igraph object
g<-graph_from_data_frame(data, vertices=vertices,directed = FALSE)
vcount(g) # number of nodes
ecount(g) # number of edges

## Set node colors by ontology category
color <- c(rgb(65,179,194,maxColorValue = 255),
           rgb(255,255,0,maxColorValue = 255),
           rgb(201,216,197,maxColorValue = 255))
names(color) <- unique(V(g)$type)
V(g)$point.col <- color[match(V(g)$type,names(color))]

## Plot GO term similarity network
plot.igraph(g,
            vertex.color=V(g)$point.col,
            shape = 1,
            vertex.size = 8,
            vertex.label.cex = 0.5,
            edge.arrow.size = 0.5)

### Community analysis: to partition the network into distinct modules
module <- cluster_edge_betweenness(g)
dendPlot(module)
plot(module,g)

## Select representative GO term for each module
module.table <- data.frame(ModuleID = membership(module)) %>%
  rownames_to_column(var = "ID") %>%
  left_join(enrich_go_ALL, by = "ID")

set.seed(2024)
PATH <- module.table %>% group_by(ModuleID) %>% slice_min(p.adjust) %>% group_by(ModuleID) %>% slice_min(pvalue) %>% group_by(ModuleID) %>% slice_max(Count) %>% group_by(ModuleID) %>% sample_n(1)

PATH$ModuleID<-as.character(PATH$ModuleID)
for (i in 1:56) {
  PATH$pathway_id[which(PATH$ModuleID==i)]<-paste0(module.table[which(module.table$ModuleID==i),c("ID")],collapse="; ")
}

## Assign remaining GO terms to "other" category
other<-enrich_go_ALL[which(!enrich_go_ALL$ID %in% module.table$ID),]
other$ModuleID<-"other_go"
other$pathway_id<-other$ID

## Combine representative and other GO terms
data<-rbind(PATH,other)
write.csv(data,file = "GO_representative_terms_protein_includ_other.csv")


#################### KEGG enrichment ############################
enrich_kegg <- enrichKEGG(gene=genes,
                          keyType = "kegg",
                          pAdjustMethod = "BH",
                          pvalueCutoff=0.05,
                          qvalueCutoff  = 0.05,
                          minGSSize = 10)

### Calculate similarity among KEGG pathways
KEGG_JC_index <- pairwise_termsim(enrich_kegg, method = "JC")@termsim # Jaccard index
sim_KEGG_corr.table <- melt(KEGG_JC_index, na.rm = TRUE)

### Select similarity score greater than 0.3 to construct similarity network
data <- sim_KEGG_corr.table %>%
  filter(value >= 0.3, Var1 != Var2)
head(data)

## Compute degree (number of connections) for each node
vertices<-c(as.character(data$Var1), as.character(data$Var2)) %>% as_tibble()
vertices<-data.frame(table(vertices$value))
colnames(vertices) <- c("node","n")

## Add KEGG pathway category
type <- enrich_kegg[,c("Description","category")]
colnames(type)<-c("node","type")

vertices %>%
  left_join(type, by="node")-> vertices

## Build igraph object
g<-graph_from_data_frame(data, vertices=vertices,directed = FALSE)
vcount(g)
ecount(g)

## Set node colors for KEGG categories
color <- c("#f57c6e","#f2b56f","#fae69e","#84c3b7","#88d8db","#71b7ed","#b8aeeb","#f2a7da","#f2a8bb")
names(color) <- unique(V(g)$type)
V(g)$point.col <- color[match(V(g)$type,names(color))]

## Plot KEGG network
plot.igraph(g,
            vertex.color=V(g)$point.col,
            shape = 1,
            vertex.size = 8,
            vertex.label.cex = 0.5,
            edge.arrow.size = 0.5)

### community analysis: to partition the network into distinct modules
module <- cluster_fast_greedy(g)
dendPlot(module)
plot(module,g)

## Select representative KEGG pathway for each module
enrich_kegg<-data.frame(enrich_kegg)
module.table <- data.frame(ModuleID = membership(module)) %>%
  rownames_to_column(var = "Description") %>%
  left_join(enrich_kegg, by = "Description")

PATH <- module.table %>% group_by(ModuleID) %>% slice_min(p.adjust) %>% group_by(ModuleID) %>% slice_min(pvalue) %>% group_by(ModuleID) %>% slice_max(Count) # %>% group_by(ModuleID) %>% sample_n(1)

PATH$ModuleID<-as.character(PATH$ModuleID)

for (i in 1:6) {
  PATH$pathway_id[which(PATH$ModuleID==i)]<-paste0(module.table[which(module.table$ModuleID==i),c("ID")],collapse="; ")
}

## Assign remaining KEGG pathways to "other" category
other<-enrich_kegg[which(!enrich_kegg$ID %in% module.table$ID),]
other$ModuleID<-"other_kegg"
other$pathway_id<-other$ID

## Combine representative and other pathways
data<-rbind(PATH,other)
write.csv(data,file = "KEGG_representative_pathways_protein_includ_other.csv")


#################### Reactome enrichment ###########################
enrich_reactome <- enrichPathway(gene= genes, 
                                 organism = "human",
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 0.05,
                                 readable = TRUE)

### Calculate similarity among Reactome pathways
reactome_JC_index <- pairwise_termsim(enrich_reactome, method = "JC")@termsim
sim_reactome_corr.table <- melt(reactome_JC_index, na.rm = TRUE)

### Select similarity score greater than 0.5 to construct similarity network
data <- sim_reactome_corr.table %>%
  filter(value >= 0.5, Var1 != Var2)
head(data)

## Compute degree (number of connections) for each node
vertices<-c(as.character(data$Var1), as.character(data$Var2)) %>% as_tibble()
vertices<-data.frame(table(vertices$value))
colnames(vertices) <- c("node","n")

## Build igraph object
g<-graph_from_data_frame(data, vertices=vertices,directed = FALSE)
vcount(g)
ecount(g)

## Plot network
plot.igraph(g,
            shape = 1,
            vertex.size = 8,
            vertex.label.cex = 0.5,
            edge.arrow.size = 0.5)

### Community analysis: to partition the network into distinct modules
module <- cluster_fast_greedy(g)
dendPlot(module)
plot(module,g)

## Select representative Reactome pathway for each module
enrich_reactome<-data.frame(enrich_reactome)
module.table <- data.frame(ModuleID = membership(module)) %>%
  rownames_to_column(var = "Description") %>%
  left_join(enrich_reactome, by = "Description")
set.seed(2024)
PATH <- module.table %>% group_by(ModuleID) %>% slice_min(p.adjust) %>% group_by(ModuleID) %>% slice_min(pvalue) %>% group_by(ModuleID) %>% slice_max(Count) %>% group_by(ModuleID) %>% sample_n(1)

PATH$ModuleID<-as.character(PATH$ModuleID)

for (i in 1:20) {
  PATH$pathway_id[which(PATH$ModuleID==i)]<-paste0(module.table[which(module.table$ModuleID==i),c("ID")],collapse="; ")
}

## Assign remaining pathways to "other" category
other<-enrich_reactome[which(!enrich_reactome$ID %in% module.table$ID),]
other$ModuleID<-"other_reactome"
other$pathway_id<-other$ID

## Combine representative and other pathways
data<-rbind(PATH,other)
write.csv(data,file = "reactome_representative_pathways_protein_includ_other.csv")


#################### Identify distinct biological functional modules ##################
PATH_go<-read.csv(file = "GO_representative_terms_protein_includ_other.csv",row.names = 1)
PATH_kegg<-read.csv(file = "KEGG_representative_pathways_protein_includ_other.csv",row.names = 1)
PATH_reactome<-read.csv(file = "reactome_representative_pathways_protein_includ_other.csv",row.names = 1)

PATH_merge<-rbind(PATH_go[,c("Description","geneID")],PATH_kegg[,c("Description","geneID")],PATH_reactome[,c("Description","geneID")])

### Convert to lists
library(purrr)
result_list <- PATH_merge %>%
  mutate(geneID_split = strsplit(geneID, "/"))

gene_list<-as.list(result_list$geneID_split)
names(gene_list)<-result_list$Description
print(result_list)

### Compute pairwise Jaccard similarity among all pathways
corr <- data.frame(Var1 = character(), Var2 = character(), Value = numeric())
for(i in 1:length(gene_list)){
  
  for(j in 1:length(gene_list)){
    
    if(i>=j){ }
    else{
      
      value = length(intersect(gene_list[[i]], gene_list[[j]]))/length(union(gene_list[[i]], gene_list[[j]]))
      
      source = names(gene_list)[i]
      
      target = names(gene_list)[j]
      
      
      corr <- rbind(corr, data.frame(Var1 = source, Var2 = target, Value = value))
      
    }
  }
}

### Select similarity score greater than 0.3 to construct the similarity network
data <- corr %>%
  filter(Value >= 0.3, Var1 != Var2)
head(data)

## node degrees
vertices<-c(as.character(data$Var1), as.character(data$Var2)) %>% as_tibble()
vertices<-data.frame(table(vertices$value))
colnames(vertices) <- c("node","n")

## Build igraph object for merged network
g<-graph_from_data_frame(data, vertices=vertices,directed = FALSE)
vcount(g) 
ecount(g)

## Plot network
plot.igraph(g,
            shape = 1,
            vertex.size = 8,
            vertex.label.cex = 0.5,
            edge.arrow.size = 0.5)

### Community analysis: to partition the network into distinct modules
module <- cluster_fast_greedy(g)
dendPlot(module)
plot(module,g)

## Assign ModuleID to each pathway
module.table <- data.frame(ModuleID = membership(module)) %>%
  rownames_to_column(var = "ID") 

## Merge p-values, adjusted p-values, counts, and gene ratios from GO, KEGG, Reactome
vertices<-PATH_merge;colnames(vertices)[1]<-"node"

# pvalue
matching_descriptions <- intersect(vertices$node, PATH_go$Description)
vertices <- vertices %>%
  left_join(
    PATH_go %>% filter(Description %in% matching_descriptions) %>% dplyr::select(Description, pvalue),
    by = c("node" = "Description")
  )

matching_descriptions <- intersect(vertices$node, PATH_kegg$Description)
vertices <- vertices %>%
  left_join(
    PATH_kegg %>% filter(Description %in% matching_descriptions) %>% dplyr::select(Description, pvalue),
    by = c("node" = "Description")
  )
vertices$pvalue <- ifelse(is.na(vertices$pvalue.x), vertices$pvalue.y, vertices$pvalue.x)
vertices <- vertices %>% dplyr::select(-c(pvalue.x,pvalue.y))

matching_descriptions <- intersect(vertices$node, PATH_reactome$Description)
vertices <- vertices %>%
  left_join(
    PATH_reactome %>% filter(Description %in% matching_descriptions) %>% dplyr::select(Description, pvalue),
    by = c("node" = "Description")
  )
vertices$pvalue <- ifelse(is.na(vertices$pvalue.x), vertices$pvalue.y, vertices$pvalue.x)
vertices <- vertices %>% dplyr::select(-c(pvalue.x,pvalue.y))

# p.adjust
matching_descriptions <- intersect(vertices$node, PATH_go$Description)
vertices <- vertices %>%
  left_join(
    PATH_go %>% filter(Description %in% matching_descriptions) %>% dplyr::select(Description, p.adjust),
    by = c("node" = "Description")
  )

matching_descriptions <- intersect(vertices$node, PATH_kegg$Description)
vertices <- vertices %>%
  left_join(
    PATH_kegg %>% filter(Description %in% matching_descriptions) %>% dplyr::select(Description, p.adjust),
    by = c("node" = "Description")
  )
vertices$p.adjust <- ifelse(is.na(vertices$p.adjust.x), vertices$p.adjust.y, vertices$p.adjust.x)
vertices <- vertices %>% dplyr::select(-c(p.adjust.x,p.adjust.y))


matching_descriptions <- intersect(vertices$node, PATH_reactome$Description)
vertices <- vertices %>%
  left_join(
    PATH_reactome %>% filter(Description %in% matching_descriptions) %>% dplyr::select(Description, p.adjust),
    by = c("node" = "Description")
  )
vertices$p.adjust <- ifelse(is.na(vertices$p.adjust.x), vertices$p.adjust.y, vertices$p.adjust.x)
vertices <- vertices %>% dplyr::select(-c(p.adjust.x,p.adjust.y))
vertices$p.adjust_log <- -log10(vertices$p.adjust)

# Count
matching_descriptions <- intersect(vertices$node, PATH_go$Description)
vertices <- vertices %>%
  left_join(
    PATH_go %>% filter(Description %in% matching_descriptions) %>% dplyr::select(Description, Count),
    by = c("node" = "Description")
  )

matching_descriptions <- intersect(vertices$node, PATH_kegg$Description)
vertices <- vertices %>%
  left_join(
    PATH_kegg %>% filter(Description %in% matching_descriptions) %>% dplyr::select(Description, Count),
    by = c("node" = "Description")
  )
vertices$Count <- ifelse(is.na(vertices$Count.x), vertices$Count.y, vertices$Count.x)
vertices <- vertices %>% dplyr::select(-c(Count.x,Count.y))

matching_descriptions <- intersect(vertices$node, PATH_reactome$Description)
vertices <- vertices %>%
  left_join(
    PATH_reactome %>% filter(Description %in% matching_descriptions) %>% dplyr::select(Description, Count),
    by = c("node" = "Description")
  )
vertices$Count <- ifelse(is.na(vertices$Count.x), vertices$Count.y, vertices$Count.x)
vertices <- vertices %>% dplyr::select(-c(Count.x,Count.y))

# gene ratio
matching_descriptions <- intersect(vertices$node, PATH_go$Description)
vertices <- vertices %>%
  left_join(
    PATH_go %>% filter(Description %in% matching_descriptions) %>% dplyr::select(Description, GeneRatio),
    by = c("node" = "Description")
  )

matching_descriptions <- intersect(vertices$node, PATH_kegg$Description)
vertices <- vertices %>%
  left_join(
    PATH_kegg %>% filter(Description %in% matching_descriptions) %>% dplyr::select(Description, GeneRatio),
    by = c("node" = "Description")
  )
vertices$GeneRatio <- ifelse(is.na(vertices$GeneRatio.x), vertices$GeneRatio.y, vertices$GeneRatio.x)
vertices <- vertices %>% dplyr::select(-c(GeneRatio.x,GeneRatio.y))

matching_descriptions <- intersect(vertices$node, PATH_reactome$Description)
vertices <- vertices %>%
  left_join(
    PATH_reactome %>% filter(Description %in% matching_descriptions) %>% dplyr::select(Description, GeneRatio),
    by = c("node" = "Description")
  )
vertices$GeneRatio <- ifelse(is.na(vertices$GeneRatio.x), vertices$GeneRatio.y, vertices$GeneRatio.x)
vertices <- vertices %>% dplyr::select(-c(GeneRatio.x,GeneRatio.y))

# pathway_id
matched_indices <- match(vertices$node, PATH_go$Description)
filtered_indices <- !is.na(matched_indices)
vertices$pathway_id[filtered_indices] <- PATH_go$pathway_id[matched_indices[filtered_indices]]

matched_indices <- match(vertices$node, PATH_kegg$Description)
filtered_indices <- !is.na(matched_indices)
vertices$pathway_id[filtered_indices] <- PATH_kegg$pathway_id[matched_indices[filtered_indices]]

matched_indices <- match(vertices$node, PATH_reactome$Description)
filtered_indices <- !is.na(matched_indices)
vertices$pathway_id[filtered_indices] <- PATH_reactome$pathway_id[matched_indices[filtered_indices]]

# ID
matched_indices <- match(vertices$node, PATH_go$Description)
filtered_indices <- !is.na(matched_indices)
vertices$ID[filtered_indices] <- PATH_go$ID[matched_indices[filtered_indices]]

matched_indices <- match(vertices$node, PATH_kegg$Description)
filtered_indices <- !is.na(matched_indices)
vertices$ID[filtered_indices] <- PATH_kegg$ID[matched_indices[filtered_indices]]

matched_indices <- match(vertices$node, PATH_reactome$Description)
filtered_indices <- !is.na(matched_indices)
vertices$ID[filtered_indices] <- PATH_reactome$ID[matched_indices[filtered_indices]]

# class
matched_indices <- match(vertices$node, PATH_go$Description)
filtered_indices <- !is.na(matched_indices)
vertices$class[filtered_indices] <- "GO"

matched_indices <- match(vertices$node, PATH_kegg$Description)
filtered_indices <- !is.na(matched_indices)
vertices$class[filtered_indices] <- "KEGG"

matched_indices <- match(vertices$node, PATH_reactome$Description)
filtered_indices <- !is.na(matched_indices)
vertices$class[filtered_indices] <- "Reactome"

## Assign remaining pathways to "other" category
other<-vertices[which(!vertices$node %in% module.table$ID),]
other$ModuleID<-"other"

## Select representative pathways for each module
vertices <- merge(vertices,module.table,by.x="node",by.y="ID")

set.seed(2024)
PATH <- vertices %>% group_by(ModuleID) %>% slice_min(p.adjust) %>% group_by(ModuleID) %>% slice_min(pvalue) %>% group_by(ModuleID) %>% slice_max(Count) %>% group_by(ModuleID) %>% sample_n(1)
PATH$ModuleID<-as.character(PATH$ModuleID)

for (i in 1:12) {
  PATH$pathway_id[which(PATH$ModuleID==i)]<-paste0(vertices[which(vertices$ModuleID==i),c("pathway_id")],collapse="; ")
  PATH$class[which(PATH$ModuleID==i)]<-paste0(vertices[which(vertices$ModuleID==i),c("class")],collapse="; ")
}

data<-rbind(PATH,other)

write.csv(data,file = "functional_module_pathways_includ_other.csv")


#################### Visualization ##################
###### Network visualization #####
### Load data
results_data<-read.csv(file = "functional_module_pathways_includ_other.csv",row.names = 1)
top_10 <- results_data[order(results_data$p.adjust), ][1:10, 1:2]

result_list <- top_10 %>%
  mutate(geneID_split = strsplit(geneID, "/"))

gene_list<-as.list(result_list$geneID_split)
names(gene_list)<-result_list$node
print(result_list)

### Compute pairwise Jaccard similarity among top 10 pathways
corr <- data.frame(Var1 = character(), Var2 = character(), Value = numeric())
for(i in 1:length(gene_list)){
  
  for(j in 1:length(gene_list)){
    
    if(i>=j){ }
    else{
      
      value = length(intersect(gene_list[[i]], gene_list[[j]]))/length(union(gene_list[[i]], gene_list[[j]]))
      
      source = names(gene_list)[i]
      
      target = names(gene_list)[j]
      
      
      corr <- rbind(corr, data.frame(Var1 = source, Var2 = target, Value = value))
      
    }
  }
}

### select similarity score greater than 0.1 to construct the similarity network
data <- corr %>%
  filter(Value >= 0.1, Var1 != Var2)
head(data)

## Compute degree (number of connections) for each node
vertices<-c(as.character(data$Var1), as.character(data$Var2)) %>% as_tibble()
vertices<-data.frame(table(vertices$value))
colnames(vertices) <- c("node","n")


## Build igraph object and plot
g<-graph_from_data_frame(data, vertices=vertices,directed = FALSE)
vcount(g)
ecount(g)

V(g)$degree <- degree(g)

plot.igraph(g,
            shape = 1,
            vertex.size = 8,
            vertex.label.cex = 0.5,
            edge.arrow.size = 0.5)

### Community analysis: to partition the network into distinct modules
module <- cluster_edge_betweenness(g)
dendPlot(module)
plot(module,g)

## Assign ModuleID to each node
module.table <- data.frame(ModuleID = membership(module)) %>%
  rownames_to_column(var = "ID") 

vertices <- merge(vertices,module.table,by.x="node",by.y="ID")
vertices <- merge(vertices,results_data[,c("node","p.adjust_log")],by="node")

### Rebuild igraph object for plotting
g<-graph_from_data_frame(data, vertices=vertices,directed = FALSE)
vcount(g)
ecount(g)

### Define node colors by ModuleID
color <- c("#478ecc",
           "#ee7e77",
           "#fdbd1a"
)
names(color) <- unique(V(g)$ModuleID)
V(g)$point.col <- color[match(V(g)$ModuleID,names(color))]

color2 <- c(rgb(199,222,248,maxColorValue = 255),
            "#e7bfc0",
            "#fdf1d1")
names(color2) <- unique(V(g)$ModuleID)
V(g)$point.col2 <- color2[match(V(g)$ModuleID,names(color2))]

### Plot
plot.igraph(g,
            vertex.color=V(g)$point.col,
            vertex.frame.color ="black",
            vertex.border=V(g)$point.col,
            shape = 1,
            vertex.size=V(g)$p.adjust_log^0.5*2, 
            rescale =TRUE,
            vertex.label=g$name,
            vertex.label.cex=1,
            vertex.label.dist=1.8,
            vertex.label.col="black",
                        mark.groups =list(V(g)$name[V(g)$ModuleID %in% names(color2)[1]],
                              V(g)$name[V(g)$ModuleID %in% names(color2)[2]],
                              V(g)$name[V(g)$ModuleID %in% names(color2)[3]]),
            mark.col=color2,
            mark.border=color2,
            edge.arrow.size=0.5,
            edge.width=abs(E(g)$Value)*10
)

###### Enrichment bubble plot ############
enrich_go<-read.csv(file = "FIsum_GO_enrich_ALL.csv")
enrich_kegg<-read.csv(file = "Enrichment/FIsum_KEGG_enrich.csv")
enrich_reactome<-read.csv(file = "Enrichment/FIsum_Reactome_enrich.csv")

top_3_go_bp <- enrich_go[enrich_go$ONTOLOGY=="BP",] %>%
  arrange(p.adjust) %>%
  select(Description, p.adjust) %>%
  mutate(dataset="GO_BP") %>%
  slice(1:3) %>%
  arrange(desc(p.adjust))
top_3_go_cc <- enrich_go[enrich_go$ONTOLOGY=="CC",] %>%
  arrange(p.adjust) %>%
  select(Description, p.adjust) %>%
  mutate(dataset="GO_CC") %>%
  slice(1:3) %>%
  arrange(desc(p.adjust))
top_3_go_mf <- enrich_go[enrich_go$ONTOLOGY=="MF",] %>%
  arrange(p.adjust) %>%
  select(Description, p.adjust) %>%
  mutate(dataset="GO_MF") %>%
  slice(1:3) %>%
  arrange(desc(p.adjust))
top_3_kegg <- enrich_kegg %>%
  arrange(p.adjust) %>%
  select(Description, p.adjust) %>%
  mutate(dataset="KEGG") %>%
  slice(1:3) %>%
  arrange(desc(p.adjust))
top_3_reactome <- enrich_reactome %>%
  arrange(p.adjust) %>%
  select(Description, p.adjust) %>%
  mutate(dataset="Reactome") %>%
  slice(1:3) %>%
  arrange(desc(p.adjust))

df<-rbind(top_3_reactome,top_3_kegg,top_3_go_mf,top_3_go_cc,top_3_go_bp)
df$p.adjust_log<- -log10(df$p.adjust)
df$dataset<-as.factor(df$dataset)
df$Description <- factor(df$Description, 
                         levels = unique(df$Description))

ggplot(df, aes(x = p.adjust_log, y = Description)) +
  geom_col(fill = "#c6c8c9", width = 0.15) +
  geom_point(aes(x = p.adjust_log, y = Description,fill = dataset,color = dataset), 
             shape = 21,
             size = 4,
             stroke = 1.5) + 
  scale_fill_manual(values = c("#527196","#478ecc","#75b5dc","#dc917b","#81b095")) +
  scale_color_manual(values = c("#527196","#478ecc","#75b5dc","#dc917b","#81b095")) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line.x = element_line(color = 'black', size = 1),
    axis.line.y = element_line(color = 'black', size = 1),
    axis.text = element_text(size = 12)
  ) +
  labs(x = '-log10(p.adjust)', y = 'pathways', title = 'Enriched Pathway') +
  theme(plot.title = element_text(hjust = 0.5, size = 14))


###### Functional module bubble plot ##############
data<-read.csv(file = "functional_module_pathways_includ_other.csv",row.names = 1)
df <- data[order(data$p.adjust), ][1:10, c("node","p.adjust","Count")]
df <- df %>% arrange(desc(p.adjust))

df$p.adjust_log<- -log10(df$p.adjust)
df$Description <- factor(df$node, 
                         levels = unique(df$node))

p = ggplot(df,aes(Count,Description))+
  geom_point(aes(size=Count,color=p.adjust_log))+
  scale_color_gradient(low="#478ecc",high = "#fc4e00")+
  labs(color=expression(-log[10](p.adjust)),size="Count",  
       x="Count",y="Pathway name",title="Pathway enrichment")+
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line.x = element_line(color = 'black', size = 1),
    axis.line.y = element_line(color = 'black', size = 1),
    axis.text = element_text(size = 12)
  )
p
