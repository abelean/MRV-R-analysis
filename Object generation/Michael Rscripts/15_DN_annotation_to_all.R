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
library( escape )

#Increase memory usage
options(future.globals.maxSize = 30000 * 1024^2) #30GB

# Set Working Directory
wd <- getwd()
# source
source( file = file.path( wd , "src" , "shared_r_functions.R" ) )

# Load the Seurat Object 'scd' ("Single Cell Data")
rds.path <- file.path( object_path , "scd.annotated.rds" )
scd <- readRDS(rds.path)
scd # Sanity check
head( Idents(object = scd) ) # Sanity check

sub.path <- file.path( object_path , "DN.scd.annotated.rds" )
dn.sub <- readRDS(sub.path)
dn.sub # Sanity check

# Create lin_3 dataframe for scd
scd.df <- data.frame( barcode = colnames(scd) , lin_1 = scd$lin_1 )
rownames( scd.df ) <- colnames( scd )
dn.df <- data.frame( barcode = colnames(dn.sub) , lin_3 = dn.sub$lin_3 )
scd.df <- left_join( scd.df , dn.df , by = 'barcode' )
for ( i in 1:nrow(scd.df) ) {
  if ( is.na(scd.df$lin_3[i]) ) {
    scd.df$lineage[i] <- scd.df$lin_1[i]
  } else {
    scd.df$lineage[i] <- scd.df$lin_3[i]
  }
}

# Add as MetaData
scd <- AddMetaData( scd, scd.df$lineage, col.name = 'lin_3')

### PRINT DN LINEAGES HIGHLIGHTED ON UMAP ###
Idents(scd) <- 'lin_3'
Idents(scd) <- factor(Idents(scd), levels =  sort( levels( Idents( scd ) ) ) ) # Needed to alphabetize the cluster level
unq <- unique( levels(Idents(scd) ) )[ ! levels( Idents(scd) ) %in% c( 'unknown' ) ]
lin.colors <- c( brewer.pal( 8 , 'Dark2') , brewer.pal( 12 , 'Paired' ) )[ 1:length(unq) ]

highlight <- list() # create a list of the lineages to highlight
for ( i in unq ) {
  highlight[[i]] <- WhichCells( scd, idents = i )
}
umap.plot <- DimPlot(scd, group.by = 'lin_3', cells.highlight = highlight, cols.highlight = lin.colors, sizes.highlight = 0.1, pt.size = 0.1) + # application of colors between DimPlot and DoHeatmap is opposite
  ggtitle( 'Main Lineages' )
file.name <- 'UMAP.lin_3.highlighted.pdf'
pdf( file = file.path(wd, 'results', 'figures',  file.name ) )
  print( umap.plot )
dev.off()

# Add escape calculations
section <- 'ORF.escape'
cat( paste0( 'Section: ' , section , '.\n' ) )
# Add MetaData of infection_cell_type
to.add <- data.frame( infection_lin_3 = paste( scd$Infection , scd$lin_3 , sep = '_' ) )
rownames( to.add ) <- colnames( scd )
scd <- AddMetaData( scd , to.add$infection_lin_3 , col.name = 'infection_lin_3' )

# set vector of gene sets 
orf.gs <- read_csv( file.path( wd , 'GS_for_escape' , 'MRV.ORF.Gene.List.csv' ) ) %>%
  pull( MRV.ORFs )

GS <- list()
GS[['MRV.ORFs']] <- orf.gs

# Set Default Assay to virus
DefaultAssay( scd ) <- 'virus'
orf.scd <- GetAssayData( scd , slot = 'counts' , assay = 'virus' )

# Calculate GSEA for scd; choose n for groups based on number of cells in smallest sample
ES <- enrichIt( obj = orf.scd , gene.sets = GS , groups = 1000 , cores = 4 ) # Groups are the number of cells to separate the enrichment calculation; cores are for parallelization
cat( 'enrich calculations complete.\n' )
# Add Enrichment scores to scd
scd <- AddMetaData( scd , ES )
cat( 'AddMetaData complete.\n' )
# Calculate Statistical Significance
mrv.sub <- subset( scd , Infection == 'MRV' )
# ES2 <- data.frame( MRV.ORFs = mrv.sub$MRV.ORFs , lin_3 = mrv.sub$lin_3 )
ES2 <- data.frame( mrv.sub[[]][ , (ncol(mrv.sub[[]]) - ncol(ES) + 1):ncol(mrv.sub[[]]) ] , lin_3 = mrv.sub$lin_3 )

output.anova <- getSignificance(ES2, group = "lin_3", fit = "ANOVA")  
x <- data.frame(Gene_Set = rownames(output.anova), output.anova)
file.name <- paste( 'ORF.sig.MRV.only.ssGSEA.ANOVA', 'tsv', sep = '.' )
write.table( x,  file.path( wd, 'results', 'Gene.Sets', file.name ), sep="\t", quote = FALSE, row.names = FALSE )
cat( 'ANOVA Significance complete.\n' )
# print violin plots of gene sets
Idents( mrv.sub ) <- 'lin_3'
Idents( mrv.sub ) <- factor( Idents( mrv.sub ) , levels =  sort( levels( Idents( mrv.sub ) ) ) ) # Needed to alphabetize the cluster level

plot_list <- VlnPlot( object = mrv.sub ,
                      features = 'MRV.ORFs' ,
                      pt.size = 0,
                      assay = 'virus' ) +
  theme( text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust = 1 ) ) 

file.name <- paste( 'ORF.sig' , 'Vln' , 'Plot' , 'pdf' , sep = '.')
pdf( file = file.path(wd, 'results', 'Gene.Sets', file.name ), width = 8, height = 11 )
  print( plot_list )
dev.off()

rm( ES , ES2 , output.anova , x )
cat( 'ANOVA VlnPlots complete.\n' )
# Calculate DEGs for each cell lineage

# NORMALIZE AND SCALL DATA ON SUBSET
DefaultAssay(scd) <- "RNA"
scd <- NormalizeData(scd) # Normalize data
all.genes <- rownames(scd)
scd <- ScaleData(scd, features = all.genes) # Scale data
cat( 'Normalize & Scale Data complete.\n' )
# Find lin_3 DEGS
section <- 'lin_3.DEGs'
Idents( scd ) <- 'Infection'
cat( paste0( 'lin_3: ' , unique( scd$lin_3 ) ) )
heatmap_list <- list()
for ( i in unique( scd$lin_3 )[! unique( scd$lin_3 ) %in% c( 'Fibroblasts' ) ] ) {
# for ( i in unique( scd$lin_3 ) ) {
  sub <- subset( scd , lin_3 == i )
  markers <- FindAllMarkers( sub , only.pos = TRUE , min.pct = 0.25 , logfc.threshold = 0.25 )
  x <- markers %>% 
#    filter( p_val_adj < 0.05 ) %>% 
    group_by(cluster) %>% 
    arrange( desc(avg_log2FC), .by_group = TRUE)
  
  file.name <- paste0( i , '.DEGs.tsv')
  write.table( x, file.path( wd, 'results', 'DEGs', file.name ), sep="\t", quote = FALSE, row.names = FALSE)
  cat( paste0( i , " DEGs saved.\n" ) )
  
  heatmap_list[[i]] <- DoHeatmap( sub , features = x$gene, group.by = 'Infection' , assay = 'SCT' ) +
    theme(axis.text.y = element_text(size = 2) ) + # decrease size of gene names
  
  rm( sub  , markers , x , file.name )
}

for ( i in unique( scd$lin_3 )[! unique( scd$lin_3 ) %in% c( 'Fibroblasts' ) ] ) {
  file.name <- paste0( i , '.heatmap.pdf')
  pdf( file = file.path(wd, 'results', 'DEGs', file.name ), width = 8, height = 11 )
    print( heatmap_list[[i]] )
  dev.off()
}
rm( heatmap_list )
cat( paste( section , ' complete.\n' ) )

# Find 'DN' in lin_1 DEGS
section <- 'lin_1.DN.DEGs'
Idents( scd ) <- 'Infection'

# subset and fin EGs for 'DN' in lin_1
section <- 'DN.DEGs'

sub <- subset( scd , lin_1 == 'DN' )
markers <- FindAllMarkers( sub , only.pos = TRUE , min.pct = 0.25 , logfc.threshold = 0.25 )
x <- filter( markers, p_val_adj < 0.05 )%>% 
  group_by(cluster) %>% 
  arrange( desc(avg_log2FC), .by_group = TRUE)
  
file.name <- paste0( 'DN' , '.DEGs.tsv')
write.table( x, file.path( wd, 'results', 'DEGs', file.name ), sep="\t", quote = FALSE, row.names = FALSE)
cat( paste0( i , " DEGs saved.\n" ) )
  
heatmap_list <- DoHeatmap( sub , features = x$gene, group.by = 'Infection' , assay = 'SCT' ) +
  theme(axis.text.y = element_text(size = 2) ) + # decrease size of gene names
  
rm( sub  , markers , x , file.name )

file.name <- paste0( 'DN' , '.heatmap.pdf')
pdf( file = file.path(wd, 'results', 'DEGs', file.name ), width = 8, height = 11 )
  print( heatmap_list )
dev.off()

cat( paste( section , ' complete.\n' ) )

# Save the object
saveRDS(scd, file = file.path( object_path , "scd.lin_3.rds" ) )
cat("Save Object complete.\n")

gc()



