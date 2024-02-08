#Make clustering and annotations for DN
#load libaries
library( ggplot2 )
library( dplyr ) 
library( tidyr )
# library( readr )
library( purrr )
library( tibble )
library( stringr )
library( Seurat )
# library( sctransform )
# library( DoubletFinder )
library( RColorBrewer )
# library( DoMultiBarHeatmap )
# library( grid )
# library( escape )
library( harmony )
#library( VGAM )
library( patchwork )
library( ggpubr )
library( enrichR )
library(reshape2)

#Increase memory usage
options(future.globals.maxSize = 32000 * 1024^2) #32GB

# source
wd <- getwd()
source( file = file.path(wd, "src", "shared_r_functions.R" ) )

# Load the Seurat Object 'tsub' ("Single Cell Data") or 'tsub' ("T cell subset")
rds.path <- file.path(wd, 'scd.lin_3.integrated.new.rds')
scd.raw <- readRDS(rds.path)
scd.raw # Sanity check
head( Idents(object = scd.raw) ) # Sanity check

DefaultAssay(scd.raw) <- "RNA"

####
scd.DN.sub <- subset(scd.raw, lin_XX %in% c('DN1','DN2','DN3','DN4'))

table(scd.DN.sub$lin_XX)
table(scd.raw$lin_XX)#Sanity check


DefaultAssay(scd.DN.sub) <- 'RNA'
scd.DN.sub[['SCT']] <- NULL
scd.DN.sub[['integrated']] <- NULL

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
scd.cc <- scd.DN.sub
scd.cc <- NormalizeData( scd.cc )
scd.cc <- FindVariableFeatures( scd.cc , selection.method = "vst" )
scd.cc <- ScaleData( scd.cc , features = rownames( scd.cc ) )
scd.cc <- CellCycleScoring(scd.cc, # min.cells = 1,
                           s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
cc.toadd <- data.frame( row.names = colnames(scd.cc) , S.Score = scd.cc$S.Score , G2M.Score = scd.cc$G2M.Score )
scd.DN.sub <- AddMetaData( object = scd.DN.sub , metadata = cc.toadd )
rm( scd.cc , cc.toadd )
cat('Cell Cycle Score complete.\n')

# Split Object and run SCTransformm while regressing out both mitochondrial genes and cell cycle genes
scd.list <- SplitObject(scd.DN.sub, split.by = "Infection") # NORMALIZE BY SAMPLE #
for (i in names(scd.list)) {
  scd.list[[i]] <- PercentageFeatureSet(scd.list[[i]], pattern = "^MT-", col.name = "percent.mt") # QC cells for mitochondiral genes
  scd.list[[i]] <- SCTransform(scd.list[[i]], vars.to.regress = c("percent.mt","S.Score", "G2M.Score"), verbose = FALSE) # run SCTransform
}
scd.features <- SelectIntegrationFeatures(object.list = scd.list, nfeatures = 3000)
cat('SelectIntegrationFeatures complete.\n')
scd.list <- PrepSCTIntegration(object.list = scd.list, anchor.features = scd.features)
cat('PrepSCTIntegration complete.\n')

# Crashed on IntegrateData (long vectors not supported yet), so switched to using a reference dataset ; source: https://satijalab.org/seurat/archive/v3.0/integration.html #
reference_dataset <- which( names(scd.list) %in% c( "Mock" ) ) # MAKE SURE THE REFERENCE DATASET IS DEFINED CORRECTLY!!! (I.E. LISTED UNDER THE PARAMETER 'SPLIT BY')
scd.anchors <- FindIntegrationAnchors(object.list = scd.list, normalization.method = "SCT",
                                      anchor.features = scd.features, reference = reference_dataset)

cat('FindIntegrationAnchors complete.\n')
scd.DN <- IntegrateData(anchorset = scd.anchors, normalization.method = "SCT", k.weight = 95)
cat('Integration complete.\n')

scd.DN <- RunPCA(scd.DN, verbose = FALSE)
scd.DN <- RunUMAP(scd.DN, reduction = "pca", dims = 1:30)

####Clustering####
#Choose clusters
scd.DN <- FindNeighbors(scd.DN, dims = 1:30, verbose = FALSE)
scd.DN <- FindClusters(scd.DN, resolution = 0.1)

# Print UMAP with chosen cluster resolution (located in meta.data under "integrated_snn_res.0.4" through "integrated_snn_res.0.6")
chosen_RES <- 'integrated_snn_res.0.1'

#Findallmarkers
Idents(scd.DN) <- "seurat_clusters"
scd.markers2 <- FindAllMarkers(scd.DN, only.pos = TRUE
                               #, min.pct = 0.25, logfc.threshold = 0.25
)

scd.markers2 <- subset(scd.markers2, p_val_adj < 0.05)
scd.markers2 %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> scd.markers2.DN

write.csv(scd.markers2.DN, file.path(wd, 'write.6', 'DN', 'DN.cluster.markers.csv'))

# Print Plots grouped by MetaData
n <- colnames(scd.DN@meta.data)[ ! colnames(scd.DN@meta.data) %in% c('RNA_CRoutput', 'TCR_CRoutput', 'BCR_CRoutput','integrated_snn_res.2.9',
                                                                     'integrated_snn_res.0.8', 'G2M.Score', 'MRV.ORFs', 'nCount_RNA',
                                                                     'nCount_SCT', 'nCount_virus', 'nFeature_RNA', 'nFeature_SCT',
                                                                     'nFeature_virus', 'S.Score', 'sex_linked') ]
plot_list <- list()
for (i in 1:length(n) ) {
  plot_list[[i]] <- DimPlot(scd.DN, group.by = n[i], pt.size = 0.025)
}
for (i in 1:length(n) ) {
  file.name <- paste('UMAP',n[i],'pdf', sep='.')
  pdf( file = file.path(wd, 'write.6', 'DN', 'MetaData', file.name ) )
  print(plot_list[[i]])
  dev.off()
}

# Read in Lineage Annotation File and Add Annotation as MetaData to Seurat Object than convert to Idents
Idents(scd.DN) <- "seurat_clusters"
lin.path <- file.path( annotation_path , "DN_annotation.tsv" )
lin_annotation <- read.delim( lin.path ) 
chosen_clusters <- data.frame( cluster = Idents(scd.DN) )
chosen_clusters$cluster <- as.numeric(as.character(chosen_clusters$cluster)) # This step was needed because the Idents were Factors that had Levels that did not correspond numerically
new.cluster.ids <- left_join( chosen_clusters, lin_annotation, by = "cluster" ) 
rownames(new.cluster.ids) <- rownames(chosen_clusters)
scd.DN <- AddMetaData( scd.DN, new.cluster.ids$lineage, col.name = 'DN')

#UMAP of new metadata
scd.DN$DN <- factor(scd.DN$DN, c("DN1/2", 'DN3', "DN4"))
colors <- c("DN1/2" = "red",
            'DN3' = "steelblue",
            "DN4" = "green")

Idents(scd.DN) <- 'DN'
plot <- UMAPPlot(scd.DN, split.by = "Infection", cols = colors)

file.name <- paste('Cluster.DN','pdf', sep='.')
pdf( file = file.path(wd, 'write.6', 'DN', file.name ) )
print( plot )
dev.off()

#Save files
saveRDS(scd.DN, file.path(wd, 'scd.DN.rds'))

gc()