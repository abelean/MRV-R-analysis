#title: "Integrate datasets"
#author: "Elisha Roberson"
#date: "2022-02-23"

# Analysis related to Tarin MRV infection data

# command-line options
args = commandArgs( trailingOnly = TRUE )
input_path = args[1]
output_path = args[2]

# libraries
library( here )
library( tidyverse )
library( Seurat )
library( glmGamPoi )

# source
source( file = here::here( "src", "shared_r_functions.R" ) )

# load raw
load( file = input_path )

## Find anchors
sc_list <- SplitObject( object = merged_raw, split.by = "Sample" )

for ( idx in 1:length( sc_list ) ) {
  if ( "virus" %in% Assays( sc_list[[ idx ]] ) ) {
    DefaultAssay( sc_list[[ idx ]] ) <- "virus"

    integrated_data <- NormalizeData( sc_list[[ idx ]],
                                      normalization.method = "CLR",
                                      margin = 2 )

    DefaultAssay( integrated_data ) <- "RNA"
  }

  sc_list[[ idx ]] <- SCTransform( sc_list[[ idx ]],
                                   method = "glmGamPoi",
                                   vars.to.regress = c( 'percent.mt', 'sex_linked' ) )
}

#data_features <- SelectIntegrationFeatures( object.list = sc_list, nfeatures = 1250 )
data_features <- SelectIntegrationFeatures( object.list = sc_list )

sc_list <- PrepSCTIntegration( object.list = sc_list,
                               anchor.features = data_features )

dataset_anchors <- FindIntegrationAnchors( object.list = sc_list,
                                           normalization.method = "SCT",
                                           anchor.features = data_features )

# integrate
integrated_data <- IntegrateData( anchorset = dataset_anchors,
                                  normalization.method = "SCT" )

save( integrated_data, file = output_path )

# Session info
Sys.time()
getwd()
sessionInfo()
