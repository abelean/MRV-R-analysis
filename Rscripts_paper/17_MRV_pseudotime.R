#Psuedotime analysis of various MRV states
#load libaries
library( ggplot2 )
library( dplyr ) 
library( tidyr )
# library( readr )
library( purrr )
library( tibble )
library( stringr )
library( Seurat )
# library( sctransform )
# library( DoubletFinder )
library( RColorBrewer )
# library( DoMultiBarHeatmap )
# library( grid )
# library( escape )
library( harmony )
library( SeuratWrappers )
library( monocle3 )
#library( VGAM )
library( patchwork )
library( ggpubr )

#Increase memory usage
options(future.globals.maxSize = 32000 * 1024^2) #32GB

# Load the Seurat Object 'tsub' ("Single Cell Data") or 'tsub' ("T cell subset")
wd <- getwd()
source( file = file.path(wd, "src", "shared_r_functions.R" ) )

rds.path <- file.path(wd, 'scd.sub.rds')
scd <- readRDS(rds.path)
scd # Sanity check
head( Idents(object = scd) ) # Sanity check

rds.path <- file.path(wd, 'Infected.scd.rds')
scd.inf <- readRDS(rds.path)
scd.inf # Sanity check
head( Idents(object = scd.inf) ) # Sanity check

#Subset MRV
scd <- subset(scd, Infection == "MRV")

#Apply labels to main object
barcodes <- data.frame("barcodes"=rownames(scd@meta.data))
infection <- data.frame("barcodes"=rownames(scd.inf@meta.data), "Infected" = scd.inf$Infected)
list.join <- left_join(barcodes, infection, by = "barcodes")

scd <- AddMetaData(scd, list.join$Infected, "Infected")

print(table(scd$Infected))
print(table(scd.inf$Infected))#Sanity check

#change none to uninfected
scd$Infected[is.na(scd$Infected)] <- 'Un-Infected'
scd$Infected[scd$Infected == 'None'] <- 'Un-Infected'
print(table(scd$Infected))#Sanity check

#Colors
set <- brewer.pal(8, "Set1")
colors <- c('DN3' = "steelblue",
            "DN4" = "green",
            "DPbla" = "purple",
            "DPre" = "orange",
            "DPsel" = "yellow",
            "CD4_SP" = "brown",
            "CD8_SP" = "pink")

####Subset Infected
scd.sub <- subset(scd, Infected == 'Infected')

# Convert Seurat Objects to CellDataSet (cds) for Monocle3

cds <- as.cell_data_set(scd.sub)
cds <- cluster_cells(cds, resolution=0.00005)#Find resolution that work best

p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)
#Subsetting partions
#integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
#cds <- as.cell_data_set(integrated.sub)
#Trajectory analysis
cds <- learn_graph(cds, use_partition = FALSE, verbose = FALSE)
plot <- plot_cells(cds,
                   color_cells_by = "cluster",
                   label_groups_by_cluster=FALSE,
                   label_leaves=TRUE,
                   label_branch_points=TRUE,
                   label_roots=TRUE,
                   label_principal_points = TRUE)#Find root here

file.name <- 'Root.Infected.pdf'
pdf(file = file.path(wd, "write.8",file.name), height = 5)
print(plot)
dev.off()

# Pseudotime
cds <- order_cells(cds, root_pr_nodes = "Y_3")#select root here, change setting later
plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "partition",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           label_branch_points=FALSE,
           trajectory_graph_color = "grey60")

#Convert back to Seurat for further analysis
integrated.sub <- as.Seurat(cds, assay = NULL)

#Graph Data (Mock vs MRV, Pseudotime, Cell Type)

Idents(integrated.sub) <- 'lin_X2'
A <- FeaturePlot(integrated.sub, "monocle3_pseudotime", split.by = 'Infected')
B <- DimPlot(integrated.sub, cols = colors) + NoLegend()
#Idents(integrated.sub) <- 'seurat_clusters'
#C <- DimPlot(integrated.sub, split.by = 'Infection')


plot <- ggarrange(A, B,
                  ncol = 1, nrow = 2)

file.name <- 'Pseudotime.Infected.pdf'
pdf(file = file.path(wd, "write.8",file.name))
print(plot)
dev.off()

####Subset Un-Infected
scd.sub <- subset(scd, Infected == 'Un-Infected')

# Convert Seurat Objects to CellDataSet (cds) for Monocle3

cds <- as.cell_data_set(scd.sub)
cds <- cluster_cells(cds, resolution=0.00005)#Find resolution that work best

p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)
#Subsetting partions
#integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
#cds <- as.cell_data_set(integrated.sub)
#Trajectory analysis
cds <- learn_graph(cds, use_partition = FALSE, verbose = FALSE)
plot <- plot_cells(cds,
                   color_cells_by = "cluster",
                   label_groups_by_cluster=FALSE,
                   label_leaves=TRUE,
                   label_branch_points=TRUE,
                   label_roots=TRUE,
                   label_principal_points = TRUE)#Find root here

file.name <- 'Root.Un-Infected.pdf'
pdf(file = file.path(wd, "write.8",file.name), height = 5)
print(plot)
dev.off()

# Pseudotime
cds <- order_cells(cds, root_pr_nodes = "Y_4")#select root here, change setting later
plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "partition",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           label_branch_points=FALSE,
           trajectory_graph_color = "grey60")

#Convert back to Seurat for further analysis
integrated.sub <- as.Seurat(cds, assay = NULL)

#Graph Data (Mock vs MRV, Pseudotime, Cell Type)

Idents(integrated.sub) <- 'lin_X2'
A <- FeaturePlot(integrated.sub, "monocle3_pseudotime", split.by = 'Infected')
B <- DimPlot(integrated.sub, cols = colors) + NoLegend()
#Idents(integrated.sub) <- 'seurat_clusters'
#C <- DimPlot(integrated.sub, split.by = 'Infection')


plot <- ggarrange(A, B,
                  ncol = 1, nrow = 2)

file.name <- 'Pseudotime.Un-Infected.pdf'
pdf(file = file.path(wd, "write.8",file.name))
print(plot)
dev.off()

####Subset Replicating
scd.sub <- subset(scd, Infected == 'Replicating')

# Convert Seurat Objects to CellDataSet (cds) for Monocle3

cds <- as.cell_data_set(scd.sub)
cds <- cluster_cells(cds, resolution=0.00005)#Find resolution that work best

p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)
#Subsetting partions
#integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
#cds <- as.cell_data_set(integrated.sub)
#Trajectory analysis
cds <- learn_graph(cds, use_partition = FALSE, verbose = FALSE)
plot <- plot_cells(cds,
                   color_cells_by = "cluster",
                   label_groups_by_cluster=FALSE,
                   label_leaves=TRUE,
                   label_branch_points=TRUE,
                   label_roots=TRUE,
                   label_principal_points = TRUE)#Find root here

file.name <- 'Root.Replicating.pdf'
pdf(file = file.path(wd, "write.8",file.name), height = 5)
print(plot)
dev.off()

# Pseudotime
cds <- order_cells(cds, root_pr_nodes = "Y_5")#select root here, change setting later
plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "partition",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           label_branch_points=FALSE,
           trajectory_graph_color = "grey60")

#Convert back to Seurat for further analysis
integrated.sub <- as.Seurat(cds, assay = NULL)

#Graph Data (Mock vs MRV, Pseudotime, Cell Type)

Idents(integrated.sub) <- 'lin_X2'
A <- FeaturePlot(integrated.sub, "monocle3_pseudotime", split.by = 'Infected')
B <- DimPlot(integrated.sub, cols = colors) + NoLegend()
#Idents(integrated.sub) <- 'seurat_clusters'
#C <- DimPlot(integrated.sub, split.by = 'Infection')


plot <- ggarrange(A, B,
                  ncol = 1, nrow = 2)

file.name <- 'Pseudotime.Replicating.pdf'
pdf(file = file.path(wd, "write.8",file.name))
print(plot)
dev.off()

gc()