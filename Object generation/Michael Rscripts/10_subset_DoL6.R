##############################################
# ANNOTATION AND SUBSETTING OF SEURAT OBJECT #
##############################################

# DoMultiBarHeatmap source: https://github.com/satijalab/seurat/issues/2201

#load libaries
library( ggplot2 )
library( dplyr ) 
library( tidyr )
library( readr )
library( purrr )
library( tibble )
library( stringr )
library( Seurat )
library( sctransform )
# library( DoubletFinder )
library( RColorBrewer )
# library( DoMultiBarHeatmap )
library( grid )


#Increase memory usage
options(future.globals.maxSize = 30000 * 1024^2) #30GB

# Set Working Directory
wd <- getwd()
# source
source( file = file.path( wd , "src" , "shared_r_functions.R" ) )

# Load the Seurat Object 'scd' ("Single Cell Data")
load('/storage1/fs1/yokoyama/Active/Tarin_MRV/2022_spring_data/results/objects/integrated_with_cluster.RData') # integrated_data is the seurat object
# Subset just Day of Life 6 ffor analysis
scd <- subset( integrated_data , Day_of_Life == 'DoL 6' )

# Print Plots of Lineage Markers #
# NEEDS A .TSV FILE OF 4 MARKERS FOR EACH LINEAGE
DefaultAssay(scd) <- "SCT"
lin.markers.path <- file.path( wd , 'Lin_markers' , 'lin_markers.tsv' )
lin.markers <- read.delim(lin.markers.path) 
n <- colnames(lin.markers)
plot_list <- list()
for (i in 1:length(n) ) {
  plot_list[[i]] <- FeaturePlot(scd, features = lin.markers[[i]], 
                                pt.size = 0.1, 
                                ncol = 2)
}

for (i in 1:length(n) ) {
  file.name <- paste( 'UMAP' , n[i] , 'pdf' , sep='.' )
  pdf( file = file.path( figure_path , 'lin_markers' , file.name ) )
    print(plot_list[[i]])
  dev.off()
}
cat("Printed UMAP based on Lineage Markers.\n")

plot_list <- list()
for (i in 1:length(n) ) {
  plot_list[[i]] <- VlnPlot(scd, features = lin.markers[[i]], 
                            pt.size = 0.1,
                            ncol = 2)
}

for (i in 1:length(n) ) {
  file.name <- paste( 'Violin' , n[i] , 'pdf' , sep='.' )
  pdf( file = file.path( figure_path , 'lin_markers', file.name ) )
    print(plot_list[[i]])
  dev.off()
}
cat("Printed Violin based on Lineage Markers.\n")

# Save the objects
rds.path <- file.path( object_path , 'scd.rds' )
saveRDS( scd , file = rds.path )

