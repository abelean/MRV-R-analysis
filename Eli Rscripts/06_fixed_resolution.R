# title: Clustering and UMAP
# author: "Elisha Roberson"
# date: "2022-03-31"

# Analysis related to Tarin MRV infection data

# libraries
library( here )
library( tidyverse )
library( cowplot )
library( Seurat )
library( gghighlight )

# source shared functions
source( file = here( "src", "shared_r_functions.R" ) )

integrated_data_path = here( 'results', 'objects', 'integrated_data.RData' )

this_resolution = fixed_resolution # imported from shared_r_functions.R

# decide graphics parameter based on OS
graphics_device_type <- case_when(
  .Platform$OS.type == "windows" ~ "windows",
  TRUE ~ "cairo"
)

image_width <- 1920
image_height <- 1024

## Cluster based on the PCs and generate UMAPs
load( file = integrated_data_path )

integrated_data <- RunPCA( integrated_data, npcs = pc_neighbor_dims )

integrated_data <- RunUMAP( integrated_data, dims = 1:pc_neighbor_dims )

integrated_data <- FindNeighbors( integrated_data, dims = 1:pc_neighbor_dims )

integrated_data <- FindClusters( integrated_data, resolution = this_resolution )

# plots
p1 <- DimPlot( integrated_data, reduction = 'umap', group.by = "Status" )
p2 <- DimPlot( integrated_data, reduction = 'umap', group.by = "seurat_clusters", label = TRUE )
p3 <- DimPlot( integrated_data, reduction = 'umap', group.by = "Infection" )

make_ggplot_jpeg( p1, here( 'results', 'figures', 'umap_by_status.jpeg' ), graphics_device_type, width = image_width, height = image_height )
make_ggplot_jpeg( p2, here( 'results', 'figures', 'umap_by_cluster.jpeg' ), graphics_device_type, width = image_width, height = image_height )
make_ggplot_jpeg( p3, here( 'results', 'figures', 'umap_by_infection.jpeg' ), graphics_device_type, width = image_width, height = image_height )
make_ggplot_jpeg( plot_grid( p1, p2 ), here( 'results', 'figures', 'umap_merged_cowplot.jpeg' ), graphics_device_type, width = image_width, height = image_height )

## Look at fraction of cells in each cluster
meta_data <- integrated_data@meta.data %>%
  select( Status, seurat_clusters )

cell_totals <- meta_data %>%
  group_by( Status ) %>%
  summarise(
    total = length( Status )
  ) %>%
  as.data.frame( . )
rownames( cell_totals ) = cell_totals$Status

cell_totals

cluster_totals <- meta_data %>%
  group_by( Status, seurat_clusters ) %>%
  summarise(
    cells = length( Status )
  )

fraction_cells <- cluster_totals %>%
  mutate( fraction = cells / cell_totals[ Status, "total" ] ) %>%
  write_tsv( x = ., file = here( 'results', "cluster_fraction_by_status.tsv" ) )

min_cluster_number <- integrated_data$seurat_clusters %>%
  as.character( . ) %>%
  as.integer( . ) %>%
  min( . )

max_cluster_number <- integrated_data$seurat_clusters %>%
  as.character( . ) %>%
  as.integer( . ) %>%
  max( . )

cluster_levels <- seq( from = min_cluster_number, to = max_cluster_number, by = 1 ) %>%
  as.character( . )

fraction_ggplot <- mutate( fraction_cells, seurat_clusters = factor( seurat_clusters, levels = cluster_levels ) ) %>%
  ggplot( data = ., mapping = aes( x = seurat_clusters, y = fraction, fill = Status ) ) +
  geom_bar( stat = "identity", color = "black", position = position_dodge() ) +
  theme_bw() +
  gg_bigger_texts +
  scale_fill_manual( values = colorBlindPalette ) +
  ylab( "Cell fraction" ) +
  xlab( "Cluster" )

make_ggplot_jpeg( fraction_ggplot, here( 'results', 'figures', 'fraction_per_cluster.jpeg' ), graphics_device_type, width = image_width, height = image_height )

# write cluster number
data.frame( resolution = this_resolution, clusters = max_cluster_number + 1 ) %>%
  write_tsv( x = ., file = here( 'results', 'cluster_number.tsv' ) )

# save cluster object
save( integrated_data, file = here( 'results', 'objects', 'integrated_with_cluster.RData' ) )

# Session info
Sys.time()
getwd()
sessionInfo()
