# title: Initial marker search
# author: "Elisha Roberson"
# date: "2022-04-05"

# Analysis related to Tarin MRV infection data

# libraries
library( here )
library( tidyverse )
library( Seurat )
# library( metap ) ## library needed for conserved markers but we don't need to load it individually

# command-line options
args = commandArgs( trailingOnly = TRUE )
this_cluster <- args[1] %>%
  as.integer( . )

# shared functions
source( file = here::here( "src", "shared_r_functions.R" ) )

# get cluster info
cluster_info <- read_tsv( file = here( 'results', 'marker_type_tests.tsv' ), col_types = c( col_integer(), col_character() ) ) %>%
  filter( Cluster == this_cluster )

test_type <- pull( cluster_info, Type )

cluster_file_string <- case_when(
  this_cluster > 9 ~ as.character( this_cluster ),
  TRUE ~ paste0( '0', as.character( this_cluster ) )
)

outfile <- paste0( "Cluster_", cluster_file_string, "_Markers.tsv" )

# Load integrated samples
load( file = here( 'results', 'objects', 'integrated_with_cluster.RData' ) )

# find markers using conserved or general test
if ( test_type == "Conserved" ) {
  FindConservedMarkers( 
    object = integrated_data,
    ident.1 = this_cluster,
	grouping.var = "Status" ) %>%
  rownames_to_column( "gene" ) %>%
  filter( max_pval < 0.10 ) %>%
  arrange( max_pval, minimump_p_val ) %>%
  write_tsv( x = ., file = here( 'results', 'lineage_markers', outfile ) )
} else if ( test_type == "General" ) {
  FindMarkers( object =  integrated_data, 
  ident.1 = this_cluster ) %>%
  rownames_to_column( "gene" ) %>%
  filter( p_val_adj < 0.10 ) %>%
  arrange( p_val_adj, p_val ) %>%
  write_tsv( x = ., file = here( 'results', 'lineage_markers', outfile ) )
} else {
  stop( "Test type not recognized" )
}

# Session info
Sys.time()

getwd()

sessionInfo()
