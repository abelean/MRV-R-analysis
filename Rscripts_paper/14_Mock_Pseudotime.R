#Pseudo time analysis of Mock
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
rds.path <- file.path(wd, 'scd.sub.rds')
scd <- readRDS(rds.path)
scd # Sanity check
head( Idents(object = scd) ) # Sanity check

#Colors
set <- brewer.pal(8, "Set1")
colors <- c('DN3' = "steelblue",
            "DN4" = "green",
            "DPbla" = "purple",
            "DPre" = "orange",
            "DPsel" = "yellow",
            "CD4_SP" = "brown",
            "CD8_SP" = "pink")

#Subset Mock
scd.sub <- subset(scd, Infection == 'Mock')

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

file.name <- 'Root.Mock.pdf'
pdf(file = file.path(wd, "write.2",file.name), height = 5)
print(plot)
dev.off()

# Pseudotime
cds <- order_cells(cds, root_pr_nodes = "Y_22")#select root here, change setting later
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
A <- FeaturePlot(integrated.sub, "monocle3_pseudotime", label.size = 0, split.by = 'Infection') + NoLegend()
B <- DimPlot(integrated.sub, cols = colors) + NoLegend()
#Idents(integrated.sub) <- 'seurat_clusters'
#C <- DimPlot(integrated.sub, split.by = 'Infection')


plot <- ggarrange(A, B,
                  ncol = 1, nrow = 2)

file.name <- 'Pseudotime.Mock.pdf'
pdf(file = file.path(wd, "write.2",file.name))
print(plot)
dev.off()

gc()