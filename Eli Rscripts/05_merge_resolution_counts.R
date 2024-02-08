# title: Merge cluster resolutions
# author: "Elisha Roberson"
# date: "2022-02-24"

# Analysis related to Tarin MRV infection data

# libraries
library( here )
library( tidyverse )

# paths
resolution_dir_path = here( 'results', 'resolution' )
output_file_path = here::here( 'results', 'clusters_per_resolution.tsv' )

list_of_dir <- dir( path = resolution_dir_path )

list_of_files <- here( 'results', 'resolution', list_of_dir, 'cluster_number.tsv' )

# load data 
cluster_count_data <- tibble()

for ( idx in 1:length( list_of_files ) ) {
	cluster_count_data <- read_tsv( file = list_of_files[ idx ]) %>%
		rbind( cluster_count_data, . )
}

write_tsv( x = cluster_count_data, file = output_file_path )

# Session info
Sys.time()
getwd()
sessionInfo()
