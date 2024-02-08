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
rds.path <- file.path( object_path , "scd.annotated.rds" )
scd <- readRDS(rds.path)
scd # Sanity check
head( Idents(object = scd) ) # Sanity check

# Subset DN cells
Idents(scd) <- 'lin_1'
DefaultAssay( scd ) <- 'RNA'
dn.sub <- subset( scd , lin_1 == 'DN' )

dn.sub[['integrated']] <- NULL # remove integrated assay
dn.sub[['SCT']] <- NULL # remove SCT assay
cat( 'Subsetting complete.\n' )

# Split Object and run SCTransformm while regressing out both mitochondrial genes and cell cycle genes 
scd.list <- SplitObject(dn.sub, split.by = "Sample") # NORMALIZE BY SAMPLE #
for (i in names(scd.list)) {
  scd.list[[i]] <- SCTransform(scd.list[[i]], vars.to.regress = c("percent.mt","sex_linked"), verbose = FALSE) # run SCTransform
}
scd.features <- SelectIntegrationFeatures(object.list = scd.list, nfeatures = 3000)
cat('SelectIntegrationFeatures complete.\n')
scd.list <- PrepSCTIntegration(object.list = scd.list, anchor.features = scd.features)
cat('PrepSCTIntegration complete.\n')

scd.anchors <- FindIntegrationAnchors(object.list = scd.list, 
                                      normalization.method = "SCT", 
                                      anchor.features = scd.features)

cat('FindIntegrationAnchors complete.\n')    
int.sub <- IntegrateData(anchorset = scd.anchors, normalization.method = "SCT")
cat('Integration complete.\n')
# Remove TCR and BCR genes from "Variable Features" in Seurat Object
VariableFeatures(int.sub) <- VariableFeatures(int.sub)[!grepl("^TRAV|^TRAJ|^TRAC|^TRBV|^TRBD|^TRBJ|^TRBC|^IGHV|^IGHD|^IGHJ|^IGHM|^IGHE|^IGHG|^IGHA|^IGLV|^IGLD|^IGLC|^IGKV|^IGKD|^IGKC", int.sub@assays$integrated@var.features)]

# Run PCA
int.sub <- RunPCA( int.sub , verbose = FALSE )
cat("RunPCA complete.\n")

# Set dims
ndims <- 30

# Perform Dimensionality Reduction with UMAP 
int.sub <- RunUMAP(int.sub, dims = 1:ndims, verbose = FALSE) 
cat("RunUMAP complete.\n")

int.sub <- FindNeighbors(int.sub, dims = 1:ndims, verbose = FALSE)
cat("FindNeighbors complete.\n")

# Find Clusters with a high cluster resolution
nclustering <- c(1:25)
for (i in nclustering) {
  RES <- (i/10)
  int.sub <- FindClusters(int.sub, resolution = RES, verbose = FALSE)
}
cat("FindClusters complete.\n")

# Print UMAP with different cluster resolution (located in meta.data under "integrated_snn_res.0.4" through "integrated_snn_res.0.6")
plot_list <- list()
for (i in nclustering) {
  RES <- (i/10)
  grp <- paste('integrated_snn_res', RES, sep = '.')
  plot_list[[i]] <- DimPlot(int.sub, group.by = grp, label = TRUE, pt.size = 0.025) + NoLegend()
}
for (i in nclustering) {
  RES <- (i/10)
  file.name <- paste( 'UMAP', 'clusters', RES, 'pdf', sep = '.')
  pdf( file = file.path(wd, 'results', 'figures', 'DN' , 'variable_clusters' , file.name ) )
    print(plot_list[[i]])
  dev.off()
}
rm( RES, file.name, grp )
cat("Printed variable clustering.\n")
# Print Plots grouped by MetaData
n <- colnames( int.sub@meta.data )
plot_list <- list()
for (i in 1:length(n) ) {
  plot_list[[i]] <- DimPlot(int.sub, group.by = n[i], pt.size = 0.025)
}
for (i in 1:length(n) ) {
  file.name <- paste('UMAP',n[i],'pdf', sep='.')
  pdf( file = file.path(wd, 'results', 'figures', 'DN' , 'UMAP' , file.name ) )
    print(plot_list[[i]])
  dev.off()
}
cat("Printed UMAP grouped by MetaData.\n")
# Print Plots of Lineage Markers #
# NEEDS A .TSV FILE OF 4 MARKERS FOR EACH LINEAGE
DefaultAssay(int.sub) <- "SCT"
lin.markers.path <- file.path( wd , 'Lin_markers' , "lin_markers.tsv" )
lin.markers <- read.delim(lin.markers.path) 
n <- colnames(lin.markers)
plot_list <- list()
for (i in 1:length(n) ) {
  plot_list[[i]] <- FeaturePlot(int.sub, features = lin.markers[[i]], 
                                pt.size = 0.1, 
                                ncol = 2)
}

for (i in 1:length(n) ) {
  file.name <- paste('UMAP',n[i],'pdf', sep='.')
  pdf( file = file.path(wd, 'results', 'figures', 'DN' , 'UMAP' , file.name ) )
    print(plot_list[[i]])
  dev.off()
}
cat("Printed UMAP based on Lineage Markers.\n")

# Save the object
file.name <- 'DN.scd.rds'
saveRDS( int.sub , file = file.path( object_path , file.name ) )

gc()



