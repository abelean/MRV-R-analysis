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
rds.path <- file.path( object_path , "DN.scd.rds" )
scd <- readRDS(rds.path)
scd # Sanity check
head( Idents(object = scd) ) # Sanity check

# Set dims
ndims <- 30

### SELECT RESOLUTION VALUE BASED ON UMAP PLOTS! ###
RES <- c( 0.7 ) ### <-- SELECT RESOLUTION VALUE BASED ON UMAP PLOTS! ###
### SELECT RESOLUTION VALUE BASED ON UMAP PLOTS! ###

# Print UMAP with chosen cluster resolution (located in meta.data under "integrated_snn_res.0.4" through "integrated_snn_res.0.6")
chosen_RES <- paste('integrated_snn_res', RES, sep = '.')
Idents( scd ) <- chosen_RES
plot_clusters <- DimPlot(scd, group.by = chosen_RES, label = TRUE, pt.size = 0.025) + NoLegend()
file.name <- paste( 'UMAP', 'chosen' , 'clusters', RES, 'pdf', sep = '.')
pdf( file = file.path(wd, 'results', 'figures' , 'DN' , file.name ) )
  print( plot_clusters )
dev.off()

rm( file.name )
cat("Printed clustering.\n")

# Print Violin Plots of Lineage Markers Based on Chosen Cluster Resolution #
# NEEDS A .TSV FILE OF 4 MARKERS FOR EACH LINEAGE
DefaultAssay(scd) <- "SCT"
lin.markers.path <- file.path( wd , 'Lin_markers' , "lin_markers.tsv" )
lin.markers <- read.delim(lin.markers.path) 
n <- colnames(lin.markers)

plot_list <- list()
for (i in 1:length(n) ) {
  plot_list[[i]] <- VlnPlot(scd, features = lin.markers[[i]], 
                            pt.size = 0.1,
                            ncol = 2)
}

for (i in 1:length(n) ) {
  file.name <- paste('Violin',n[i],'pdf', sep='.')
  pdf( file = file.path(wd, 'results', 'figures' , 'DN' , 'Violin' , file.name ), width = 20, height = 20)
  print(plot_list[[i]])
  dev.off()
}
cat("Printed Violin based on Lineage Markers.\n")

# Print highlighted clusters
Idents(scd) <- chosen_RES
plot_list <- list()
for (i in levels( Idents( scd ) ) ) {
  highlight <- WhichCells( scd, idents = i )
  plot_list[[i]] <- DimPlot(scd, group.by = 'ident', cells.highlight = highlight, cols.highlight = '#0000FF', sizes.highlight = 0.1, pt.size = 0.025) +
    ggtitle( paste0('Cluster ', i ) ) +
    theme(legend.position = 'none')
}
for (i in levels( Idents( scd ) ) ) {
  file.name <- paste('Cluster', i, 'pdf', sep='.')
  pdf( file = file.path(wd, 'results', 'figures' , 'DN' , 'chosen_clusters', file.name ) )
  print(plot_list[[i]])
  dev.off()
}
rm( plot_list, highlight )

# Annotate DN Clusters #

# Read in Lineage Annotation File and Add Annotation as MetaData to Seurat Object than convert to Idents
lin.path <- file.path( annotation_path , "lin_annotation_3.tsv" )
lin_annotation <- read.delim( lin.path ) 
chosen_clusters <- data.frame( cluster = Idents(scd) )
chosen_clusters$cluster <- as.numeric(as.character(chosen_clusters$cluster)) # This step was needed because the Idents were Factors that had Levels that did not correspond numerically
new.cluster.ids <- left_join( chosen_clusters, lin_annotation, by = "cluster" ) 
rownames(new.cluster.ids) <- rownames(chosen_clusters)
scd <- AddMetaData( scd, new.cluster.ids$lineage, col.name = 'lin_3')

### PRINT DN LINEAGES HIGHLIGHTED ON UMAP ###
Idents(scd) <- 'lin_3'
Idents(scd) <- factor(Idents(scd), levels =  sort( levels( Idents( scd ) ) ) ) # Needed to alphabetize the cluster level
unq <- unique( levels(Idents(scd) ) )
lin.colors <- brewer.pal( 4 , 'Set1')

highlight <- list() # create a list of the lineages to highlight
for ( i in unq ) {
  highlight[[i]] <- WhichCells( scd, idents = i )
}
umap.plot <- DimPlot(scd, group.by = 'lin_3', cells.highlight = highlight, cols.highlight = lin.colors, sizes.highlight = 0.1, pt.size = 0.1) + # application of colors between DimPlot and DoHeatmap is opposite
  ggtitle( 'DN Lineages' )
file.name <- 'UMAP.lin_3.highlighted.pdf'
pdf( file = file.path(wd, 'results', 'figures', 'DN',  file.name ) )
  print( umap.plot )
dev.off()

# Save the object
saveRDS(scd, file = file.path( object_path , "DN.scd.annotated.rds" ) )
cat("Save Object complete.\n")


gc()



