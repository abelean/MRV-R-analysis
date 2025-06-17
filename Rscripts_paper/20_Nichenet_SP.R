#Nichenet analysis using SP cells as recievers and mTECs/DCs as senders, Mock vs. MRV
#load libaries
library( ggplot2 )
library( cowplot )
library( dplyr ) 
library( tidyr )
library( purrr )
library( tibble )
library( stringr )
library( Seurat )
library( pheatmap )
library( RColorBrewer )
library( VGAM )
library( readr )
library( nichenetr )
library( future )
library( readxl )

#Increase memory usage and set plan
options(future.globals.maxSize = 32 * 1024^2) #32GB
#plan("multicore")

#ID project directory
wd <- getwd()
source( file = file.path(wd, "src", "shared_r_functions.R" ) )

# Load the Seurat Object 'scd' ("Single Cell Data")
rds.path <- file.path(wd, 'scd.lin_3.integrated.new.rds')
scd <- readRDS(rds.path)
scd # Sanity check
head( Idents(object = scd) ) # Sanity check

##Assign name
name <- 'SP'

scd <- subset(scd, lin_X3 %in% c('CD4_SP','CD8_SP','mTEC','DC'))

##Combine CD4 and CD8
for(i in 1:length(scd$lin_X3)){
  if(scd$lin_X3[i] %in% c('CD4_SP','CD8_SP')){
    scd$lin_X3[i] <- "SP"
  }
}

scd$lin_X3 <- factor(scd$lin_X3, c("SP", "mTEC", "DC"))
Idents(scd) <- "lin_X3"
print(table(scd$lin_X3))#Sanity check

####Subset by Mock####
scd.condition <- subset(scd, Infection == "Mock")

###Nichenet
rds.path <- file.path( object_path , 'ligand_target_matrix.rds')
ligand_target_matrix <- readRDS(rds.path)
rds.path <- file.path( object_path , 'lr_network.rds')
lr_network <- readRDS(rds.path)
rds.path <- file.path( object_path , 'weighted_networks.rds')
weighted_networks <- readRDS(rds.path)
weighted_networks_lr = weighted_networks$lr_sig %>%
  inner_join(lr_network %>%
               distinct(from,to), by = c("from","to"))

##convert to mouse
#matrix
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols() 
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols() 

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

#lr_network
lr_network$from = lr_network$from %>% convert_human_to_mouse_symbols() 
lr_network$to= lr_network$to %>% convert_human_to_mouse_symbols()

lr_network <- na.omit(lr_network)

#weight
weighted_networks_lr$from = weighted_networks_lr$from %>% convert_human_to_mouse_symbols() 
weighted_networks_lr$to= weighted_networks_lr$to %>% convert_human_to_mouse_symbols()

weighted_networks_lr <- na.omit(weighted_networks_lr)

#Define receiver and sender
Idents(scd.condition) <- 'lin_X3'
receiver = "SP"
expressed_genes_receiver = get_expressed_genes(receiver, scd.condition, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>%
  .[. %in% rownames(ligand_target_matrix)]

sender_celltypes = c('mTEC','DC')
list_expressed_genes_sender = sender_celltypes %>%
  unique() %>%
  lapply(get_expressed_genes, scd.condition, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>%
  unlist() %>%
  unique()

#Define gene set
seurat_obj_receiver= subset(scd, idents = receiver)
Idents(seurat_obj_receiver) <- "Infection"

condition_oi = "Mock"
condition_reference = "MRV"

DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>%
  
  rownames_to_column("gene")

geneset_oi = DE_table_receiver %>%
  filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>%
  pull(gene)
geneset_oi = geneset_oi %>%
  .[. %in% rownames(ligand_target_matrix)]

#Define potential ligands
ligands = lr_network %>%
  pull(from) %>%
  unique()
receptors = lr_network %>%
  pull(to) %>%
  unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>%
  filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
  pull(from) %>%
  unique()

#Nichenet ligand analysis
ligand_activities = predict_ligand_activities(geneset = geneset_oi,
                                              background_expressed_genes = background_expressed_genes,
                                              ligand_target_matrix = ligand_target_matrix,
                                              potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>%
  arrange(-pearson) %>%
  mutate(rank = rank(desc(pearson)))
ligand_activities

best_upstream_ligands = ligand_activities %>%
  top_n(20, pearson) %>%
  arrange(-pearson) %>%
  pull(test_ligand) %>%
  unique() %>%
  str_sort(numeric = TRUE)

as.character(best_upstream_ligands) %in% rownames(scd.condition@assays$RNA@counts)

#Infer receptors
active_ligand_target_links_df = best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>%
  bind_rows() %>%
  drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df,
                                                                 ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>%
  rev() %>%
  make.names()

order_targets = active_ligand_target_links_df$target %>%
  unique() %>%
  intersect(rownames(active_ligand_target_links)) %>%
  make.names()

rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>%
  make.names() # make.names() for heatmap visualization of genes like H2-T23

colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>%
  make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>%
  t()

p_ligand_target_network = vis_ligand_target %>%
  
  make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",
                      legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + 
  theme(axis.text.x = element_text(face = "italic")) + 
  scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))

p_ligand_target_network

#receptors of top-ranked ligands
lr_network_top = lr_network %>%
  filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>%
  distinct(from,to)

best_upstream_receptors = lr_network_top %>%
  pull(to) %>%
  unique()

lr_network_top_df_large = weighted_networks_lr %>%
  filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>%
  spread("from","weight",fill = 0)

lr_network_top_matrix = lr_network_top_df %>%
  select(-to) %>%
  as.matrix() %>%
  magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>%
                      t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>%
  intersect(rownames(lr_network_top_matrix))

order_ligands_receptor = order_ligands_receptor %>%
  intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>%
  make.names()

colnames(vis_ligand_receptor_network) = order_ligands_receptor %>%
  make.names()

p_ligand_receptor_network = vis_ligand_receptor_network %>%
  t() %>%
  make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred",
                      x_axis_position = "top",legend_title = "Prior interaction potential")

p_ligand_receptor_network

lr_network_strict = lr_network %>%
  filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>%
  pull(from) %>% unique()

receptors_bona_fide = lr_network_strict %>%
  pull(to) %>%
  unique()

lr_network_top_df_large_strict = lr_network_top_df_large %>%
  distinct(from,to) %>%
  inner_join(lr_network_strict, by = c("from","to")) %>%
  distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>%
  inner_join(lr_network_top_df_large, by = c("from","to"))

lr_network_top_df_strict = lr_network_top_df_large_strict %>%
  spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>%
  select(-to) %>%
  as.matrix() %>%
  magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>%
                      t(), method = "binary")

hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>%
  intersect(rownames(lr_network_top_matrix_strict))

order_ligands_receptor = order_ligands_receptor %>%
  intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>%
  make.names()

colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>%
  make.names()

p_ligand_receptor_network_strict = lr_network_top_matrix_strict %>%
  t() %>%
  make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred",
                      x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)") +
  scale_y_discrete(limits=rev) +
  theme(text=element_text(size=20), legend.position = "right") +
  scale_fill_gradient(low = "white", high =  "mediumvioletred")
p_ligand_receptor_network_strict

file.name <- paste(name,'Mock.3','pdf', sep='.')
pdf( file = file.path( wd , 'write.9', file.name ) , width = 10 , height = 10 )
print(p_ligand_receptor_network_strict)
dev.off()

scd.condition@active.ident<- factor(scd.condition@active.ident, levels = c('SP','mTEC','DC'))##Add levels here
scd.condition <- subset(scd.condition, lin_X3 != receiver)
plot <- DotPlot(scd.condition, features = rev(colnames(lr_network_top_matrix_strict)), cols = "RdYlBu",
                scale.min = 0, scale.max = 100) + coord_flip() +
  theme(axis.text = element_text(size = 15,angle = 45, hjust = 1), legend.position = "right")

file.name <- paste(name,'Mock.0','pdf', sep='.')
pdf( file = file.path( wd , 'write.9', file.name ) , width = 5 , height = 10 )
print(plot)
dev.off()

####Subset by MRV####
scd.condition <- subset(scd, Infection == "MRV")

###Nichenet
rds.path <- file.path( object_path , 'ligand_target_matrix.rds')
ligand_target_matrix <- readRDS(rds.path)
rds.path <- file.path( object_path , 'lr_network.rds')
lr_network <- readRDS(rds.path)
rds.path <- file.path( object_path , 'weighted_networks.rds')
weighted_networks <- readRDS(rds.path)
weighted_networks_lr = weighted_networks$lr_sig %>%
  inner_join(lr_network %>%
               distinct(from,to), by = c("from","to"))

##convert to mouse
#matrix
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols() 
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols() 

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

#lr_network
lr_network$from = lr_network$from %>% convert_human_to_mouse_symbols() 
lr_network$to= lr_network$to %>% convert_human_to_mouse_symbols()

lr_network <- na.omit(lr_network)

#weight
weighted_networks_lr$from = weighted_networks_lr$from %>% convert_human_to_mouse_symbols() 
weighted_networks_lr$to= weighted_networks_lr$to %>% convert_human_to_mouse_symbols()

weighted_networks_lr <- na.omit(weighted_networks_lr)

#Define receiver and sender
Idents(scd.condition) <- 'lin_X3'
receiver = "SP"
expressed_genes_receiver = get_expressed_genes(receiver, scd.condition, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>%
  .[. %in% rownames(ligand_target_matrix)]

sender_celltypes = c('mTEC','DC')
list_expressed_genes_sender = sender_celltypes %>%
  unique() %>%
  lapply(get_expressed_genes, scd.condition, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>%
  unlist() %>%
  unique()

#Define gene set
seurat_obj_receiver= subset(scd, idents = receiver)
Idents(seurat_obj_receiver) <- "Infection"

condition_oi = "MRV"
condition_reference = "Mock"

DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>%
  
  rownames_to_column("gene")

geneset_oi = DE_table_receiver %>%
  filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>%
  pull(gene)
geneset_oi = geneset_oi %>%
  .[. %in% rownames(ligand_target_matrix)]

#Define potential ligands
ligands = lr_network %>%
  pull(from) %>%
  unique()
receptors = lr_network %>%
  pull(to) %>%
  unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>%
  filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
  pull(from) %>%
  unique()

#Nichenet ligand analysis
ligand_activities = predict_ligand_activities(geneset = geneset_oi,
                                              background_expressed_genes = background_expressed_genes,
                                              ligand_target_matrix = ligand_target_matrix,
                                              potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>%
  arrange(-pearson) %>%
  mutate(rank = rank(desc(pearson)))
ligand_activities

best_upstream_ligands = ligand_activities %>%
  top_n(20, pearson) %>%
  arrange(-pearson) %>%
  pull(test_ligand) %>%
  unique() %>%
  str_sort(numeric = TRUE)

as.character(best_upstream_ligands) %in% rownames(scd.condition@assays$RNA@counts)

#Infer receptors
active_ligand_target_links_df = best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>%
  bind_rows() %>%
  drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df,
                                                                 ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>%
  rev() %>%
  make.names()

order_targets = active_ligand_target_links_df$target %>%
  unique() %>%
  intersect(rownames(active_ligand_target_links)) %>%
  make.names()

rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>%
  make.names() # make.names() for heatmap visualization of genes like H2-T23

colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>%
  make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>%
  t()

p_ligand_target_network = vis_ligand_target %>%
  
  make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",
                      legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + 
  theme(axis.text.x = element_text(face = "italic")) + 
  scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))

p_ligand_target_network

#receptors of top-ranked ligands
lr_network_top = lr_network %>%
  filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>%
  distinct(from,to)

best_upstream_receptors = lr_network_top %>%
  pull(to) %>%
  unique()

lr_network_top_df_large = weighted_networks_lr %>%
  filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>%
  spread("from","weight",fill = 0)

lr_network_top_matrix = lr_network_top_df %>%
  select(-to) %>%
  as.matrix() %>%
  magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>%
                      t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>%
  intersect(rownames(lr_network_top_matrix))

order_ligands_receptor = order_ligands_receptor %>%
  intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>%
  make.names()

colnames(vis_ligand_receptor_network) = order_ligands_receptor %>%
  make.names()

p_ligand_receptor_network = vis_ligand_receptor_network %>%
  t() %>%
  make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred",
                      x_axis_position = "top",legend_title = "Prior interaction potential")

p_ligand_receptor_network

lr_network_strict = lr_network %>%
  filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>%
  pull(from) %>% unique()

receptors_bona_fide = lr_network_strict %>%
  pull(to) %>%
  unique()

lr_network_top_df_large_strict = lr_network_top_df_large %>%
  distinct(from,to) %>%
  inner_join(lr_network_strict, by = c("from","to")) %>%
  distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>%
  inner_join(lr_network_top_df_large, by = c("from","to"))

lr_network_top_df_strict = lr_network_top_df_large_strict %>%
  spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>%
  select(-to) %>%
  as.matrix() %>%
  magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>%
                      t(), method = "binary")

hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>%
  intersect(rownames(lr_network_top_matrix_strict))

order_ligands_receptor = order_ligands_receptor %>%
  intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>%
  make.names()

colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>%
  make.names()

p_ligand_receptor_network_strict = lr_network_top_matrix_strict %>%
  t() %>%
  make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred",
                      x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)") +
  scale_y_discrete(limits=rev) +
  theme(text=element_text(size=20), legend.position = "right") +
  scale_fill_gradient(low = "white", high =  "mediumvioletred")
p_ligand_receptor_network_strict

file.name <- paste(name,'MRV.3','pdf', sep='.')
pdf( file = file.path( wd , 'write.9', file.name ) , width = 10 , height = 10 )
print(p_ligand_receptor_network_strict)
dev.off()

scd.condition@active.ident <- factor(scd.condition@active.ident, levels = c('SP','mTEC','DC'))##Add levels here
scd.condition <- subset(scd.condition, lin_X3 != receiver)
plot <- DotPlot(scd.condition, features = rev(colnames(lr_network_top_matrix_strict)), cols = "RdYlBu",
                scale.min = 0, scale.max = 100) + coord_flip() +
  theme(axis.text = element_text(size = 15,angle = 45, hjust = 1), legend.position = "right")

file.name <- paste(name,'MRV.0','pdf', sep='.')
pdf( file = file.path( wd , 'write.9', file.name ) , width = 5 , height = 10 )
print(plot)
dev.off()


gc()
