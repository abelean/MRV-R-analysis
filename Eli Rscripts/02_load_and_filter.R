# title: "Load and filter cells"
# author: "Elisha Roberson"
# date: "2022-02-22"

# Analysis related to Tarin MRV infection data

# command-line options
args = commandArgs( trailingOnly = TRUE )
merged_output_path <- args[1]
cell_count_path <- args[2]

# load libraries
library( here )
library( tidyverse )
library( Seurat )

# Source
source( file = here::here( "src", "shared_r_functions.R" ) )

# project info
project_info <- read_tsv( file = here::here( "info", "project_info.txt" ) )

# Y chromosome gene import
all_y_chrom_genes <- read_tsv( file = here::here( 'info', 'mouse_y_genes.txt.gz' ) ) %>%
  pull( symbol ) %>%
  unique( . )

# viral ORFs
viral_orf_list <- read_tsv( file = here::here( 'info', 'viral_orfs.txt.gz' ) ) %>%
  pull( gene_id ) %>%
  unique( . )

# Read raw data
# Filter features / mito
# Remove doublets

seurat_object_name_list = character( 0 )
seurat_cell_id_list = character( 0 )

cell_numbers = data.frame( rgsm = project_info$rgsm, raw_cells = 0, filtered_cells = 0 )
rownames( cell_numbers ) = cell_numbers$rgsm

for ( row_idx in 1:nrow( project_info ) ) {
  # current sample info
  folder_name <- project_info %>%
    pull( data_folder ) %>%
    .[ row_idx ]

  rgid <- project_info %>%
    pull( rgid ) %>%
    .[ row_idx ]

  rgsm <- project_info %>%
    pull( rgsm ) %>%
    .[ row_idx ]

  infection <- project_info %>%
    pull( infection ) %>%
    .[ row_idx ]

  status <- project_info %>%
    pull( status ) %>%
    .[ row_idx ]

  day_of_life <- project_info %>%
    pull( day_of_life ) %>%
    .[ row_idx ]

  doublet_filename <- paste0( rgid, "_doublets.tsv.gz" )

  # load doublet list
  doublet_cells <- read_tsv( file = here::here( "results", "doublets", doublet_filename ) ) %>%
    pull( cell_id )

  # load data
  count_data <- Read10X( data.dir = here::here( "data", folder_name ), strip.suffix = TRUE )

  count_meta <- data.frame( cell_name = colnames( count_data ),
                            RGID = rgid,
                            Status = status,
                            Sample = rgsm,
                            Infection = infection,
                            Day_of_Life = day_of_life ) %>%
  column_to_rownames( "cell_name" )

  # convert to seurat
  seurat_data <- CreateSeuratObject( counts = count_data,
                                     project = project_name,
                                     min.cells = min_cell_cutoff,
                                     min.features = min_feature_cutoff,
                                     meta.data = count_meta )
  # raw counts
  cell_numbers[ rgsm, "raw_cells" ] = length( colnames( seurat_data ) )

  # remove intermediate objects
  rm( count_data )
  rm( count_meta )

  # Basic feature / mitochondrial QC
  gene_name_list <- rownames( seurat_data )

  # get search string for mito genes and get percentages
  mt_string <- guess_mitochondrial_gene_pattern( gene_name_list )

  if ( is.na( mt_string ) ) {
    stop( paste0( "No mitochondrial string detected" ) )
  }

  seurat_data[[ "percent.mt" ]] <- PercentageFeatureSet( seurat_data, pattern = mt_string )

  # get list of genes for the Y chromosome
  y_features_in_this_sample <- intersect( x = gene_name_list, y = all_y_chrom_genes )

  if ( length( y_features_in_this_sample ) == 0 ) {
    stop( "No Y chromosome genes found!!!" )
  }

  seurat_data[[ "sex_linked" ]] <- PercentageFeatureSet( seurat_data, features = y_features_in_this_sample )

  # filter based on mito counts
  mt_cutoff_value <- case_when(
    length( which( seurat_data[[ "percent.mt" ]][,1] > 1.0 ) ) > 0 ~ 20.0,
    TRUE ~ 0.20
  )

  seurat_data <- subset( seurat_data, nFeature_RNA > min_feature_cutoff & percent.mt < mt_cutoff_value )

  # Remove doublets
  seurat_data <- seurat_data[ , setdiff( colnames( seurat_data ), doublet_cells ) ]

  # count filtered
  cell_numbers[ rgsm, "filtered_cells" ] = length( colnames( seurat_data ) )

  # add viral transcripts if they exist
  if ( rgid %in% names( has_viral_tags ) ) {
    vdirectory <- here( 'data', has_viral_tags[ rgid ] ) %>%
      unname( . )

    cells_to_keep <- colnames( seurat_data )

    count_data <- Read10X( data.dir = here( vdirectory ), strip.suffix = TRUE )

    keep_features <- intersect( x = rownames( count_data ), y = viral_orf_list )

    count_data <- count_data[ keep_features, cells_to_keep ]

    viral_seurat_data <- CreateAssayObject( counts = count_data )

    seurat_data[[ 'virus' ]] = viral_seurat_data

    DefaultAssay( seurat_data ) <- "RNA"
  }

  # save objects and move on
  object_name <- paste0( 'seurat_', row_idx )

  seurat_object_name_list[ row_idx ] = object_name
  seurat_cell_id_list[ row_idx ] = rgsm

  assign( x = object_name, value = seurat_data )
  rm( seurat_data )
}

# merge into seurat
seurat_object_list <- list()
for ( idx in 2:length( seurat_object_name_list ) ) {
  seurat_object_list[ idx - 1 ] <- base::get( seurat_object_name_list[ idx ] )
}

merged_raw <- merge( x = base::get( seurat_object_name_list[1] ),
                     y = seurat_object_list,
                     add.cell.ids = seurat_cell_id_list,
                     project = project_name )

# save R data objects
save( merged_raw, file = merged_output_path )

# write cell numbers
write_tsv( x = cell_numbers, file = cell_count_path )

# Session info
Sys.time()
getwd()
sessionInfo()
