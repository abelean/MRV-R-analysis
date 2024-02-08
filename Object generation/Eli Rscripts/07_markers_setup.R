# title: Initial marker setup
# author: "Elisha Roberson"
# date: "2022-04-05"

# Analysis related to Tarin MRV infection data

# libraries
library( here )
library( tidyverse )
library( reshape2 )
library( Seurat )

# shared functions
source( file = here( "src", "shared_r_functions.R" ) )

# Load integrated samples
load( file = here( 'results', 'objects', 'integrated_with_cluster.RData' ) )

# cluster info & decide how to find markers
min_cluster_number <- integrated_data$seurat_clusters %>%
  as.character( . ) %>%
  as.integer( . ) %>%
  min( . )

max_cluster_number <- integrated_data$seurat_clusters %>%
  as.character( . ) %>%
  as.integer( . ) %>%
  max( . )

data.frame( Name = c( 'Min', 'Max' ), Value = c( min_cluster_number, max_cluster_number ) ) %>%
  write_tsv( x = ., file = here( 'results', 'min_max_clusters.tsv' ) )

cluster_set <- seq( min_cluster_number, max_cluster_number, by = 1 )

too_few_cell_cluster <- read_tsv( file = here( "results", "cluster_fraction_by_status.tsv" ), col_types = c( col_character(), col_integer(), col_integer(), col_double() ) ) %>%
  select( -fraction ) %>%
  arrange( seurat_clusters, Status ) %>%
  mutate( seurat_clusters = paste0( 'cluster_', seurat_clusters ) ) %>%
  mutate( seurat_clusters = fct_inorder( seurat_clusters ) ) %>%
  pivot_wider( names_from = Status, values_from = cells, values_fill = list( "cells" = 0 ) ) %>%
  melt( . ) %>%
  filter( value < 3 ) %>%
  mutate( seurat_clusters = as.character( seurat_clusters ) ) %>%
  mutate( seurat_clusters = str_replace( string = seurat_clusters, pattern = "cluster_", replacement = "" ) ) %>%
  pull( seurat_clusters ) %>%
  as.integer( . )
  
marker_type_tests <- tibble()

if ( length( too_few_cell_cluster ) == 0 ) {
  marker_type_tests <- tibble( Cluster = cluster_set, Type = "Conserved" )
} else if ( length( too_few_cell_cluster ) == length( cluster_set ) ) {
  marker_type_tests <- tibble( Cluster = cluster_set, Type = "General" )
} else {
  marker_type_tests <- setdiff( x = cluster_set, y = too_few_cell_cluster ) %>%
    tibble( Cluster = ., Type = "Conserved" )
  
  marker_type_tests <- tibble( Cluster = too_few_cell_cluster, Type = "General" ) %>%
    rbind( marker_type_tests, . ) %>%
    arrange( Cluster )
}

write_tsv( x = marker_type_tests, file = here( "results", "marker_type_tests.tsv" ) )

# Session info
Sys.time()

getwd()

sessionInfo()
