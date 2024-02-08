##############################################
# ANNOTATION AND SUBSETTING OF SEURAT OBJECT #
##############################################

#load libaries
library( ggplot2 )
library( dplyr ) 
library( tidyr )
# library( readr )
# library( purrr )
library( tibble )
# library( stringr )
library( Seurat )
library( sctransform )
# library( DoubletFinder )

#Increase memory usage
options(future.globals.maxSize = 60000 * 1024^2) #60GB

# Set Working Directory
wd <- getwd()
# source
source( file = file.path(wd, "src", "shared_r_functions.R" ) )

# Load the Seurat Object 'scd' ("Single Cell Data")
rds.path <- file.path( object_path , "scd.rds" )
scd <- readRDS(rds.path)
scd # Sanity check
head( Idents(object = scd) ) # Sanity check

# Read in Lineage Annotation File and Add Annotation as MetaData to Seurat Object than convert to Idents
lin.path <- file.path( annotation_path , "lin_annotation_1.tsv" )
lin_annotation <- read.delim( lin.path ) 
chosen_clusters <- data.frame( cluster = Idents(scd) )
chosen_clusters$cluster <- as.numeric(as.character(chosen_clusters$cluster)) # This step was needed because the Idents were Factors that had Levels that did not correspond numerically
new.cluster.ids <- left_join( chosen_clusters, lin_annotation, by = "cluster" ) 
rownames(new.cluster.ids) <- rownames(chosen_clusters)
scd <- AddMetaData( scd, new.cluster.ids$lineage, col.name = 'lin_1')

rm( lin.path , lin_annotation , new.cluster.ids )

# Repeat for Lineage_2
lin.path <- file.path( annotation_path , "lin_annotation_2.tsv" )
lin_annotation <- read.delim( lin.path ) 
new.cluster.ids <- left_join( chosen_clusters, lin_annotation, by = "cluster" ) 
rownames(new.cluster.ids) <- rownames(chosen_clusters)
scd <- AddMetaData( scd, new.cluster.ids$lineage, col.name = 'lin_2')

### PRINT MAIN LINEAGES HIGHLIGHTED ON UMAP ###
Idents(scd) <- 'lin_1'
Idents(scd) <- factor(Idents(scd), levels =  sort( levels( Idents( scd ) ) ) ) # Needed to alphabetize the cluster level
unq <- unique( levels(Idents(scd) ) )[ ! levels( Idents(scd) ) %in% c( 'unknown' ) ]
lin.colors <- c( brewer.pal( 8 , 'Dark2') , brewer.pal( 12 , 'Paired' ) )[ 1:length(unq) ]

highlight <- list() # create a list of the lineages to highlight
for ( i in unq ) {
    highlight[[i]] <- WhichCells( scd, idents = i )
}
umap.plot <- DimPlot(scd, group.by = 'lin_1', cells.highlight = highlight, cols.highlight = lin.colors, sizes.highlight = 0.1, pt.size = 0.025) + # application of colors between DimPlot and DoHeatmap is opposite
        ggtitle( 'Main Lineages' )
file.name <- 'UMAP.lin_1.highlighted.pdf'
pdf( file = file.path(wd, 'results', 'figures',  file.name ) )
    print( umap.plot )
dev.off()
rm( highlight )

# Repeat for lin_2
Idents(scd) <- 'lin_2'
Idents(scd) <- factor(Idents(scd), levels =  sort( levels( Idents( scd ) ) ) ) # Needed to alphabetize the cluster level
unq <- unique( levels(Idents(scd) ) )[ ! levels( Idents(scd) ) %in% c( 'unknown' ) ]
lin.colors <- c( brewer.pal( 8 , 'Dark2') , brewer.pal( 12 , 'Paired' ) )[ 1:length(unq) ]

highlight <- list() # create a list of the lineages to highlight
for ( i in unq ) {
  highlight[[i]] <- WhichCells( scd, idents = i )
}
umap.plot <- DimPlot(scd, group.by = 'lin_2', cells.highlight = highlight, cols.highlight = lin.colors, sizes.highlight = 0.1, pt.size = 0.025) + # application of colors between DimPlot and DoHeatmap is opposite
  ggtitle( 'Main Lineages' )
file.name <- 'UMAP.lin_2.highlighted.pdf'
pdf( file = file.path(wd, 'results', 'figures',  file.name ) )
print( umap.plot )
dev.off()
rm( highlight )
### END OF SECTION ###

# Save the object
saveRDS( scd , file = file.path( object_path , "scd.annotated.rds" ) )

gc()
