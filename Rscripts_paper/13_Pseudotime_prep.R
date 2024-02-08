#Generating dimenstionality reduction of main clusters for pseudo time
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
library( SeuratWrappers )
library( monocle3 )
#library( VGAM )
library( patchwork )
library( ggpubr )

#Increase memory usage
options(future.globals.maxSize = 32000 * 1024^2) #32GB

# Load the Seurat Object 'tsub' ("Single Cell Data") or 'tsub' ("T cell subset")
wd <- getwd()
rds.path <- file.path(wd, 'scd.lin_3.integrated.new.rds')
tsub <- readRDS(rds.path)
tsub # Sanity check
head( Idents(object = tsub) ) # Sanity check

#subset out smaller clusters
scd.raw <- subset(tsub, lin_X2 %in% c("DN3","DN4","DPsel","DPre",
                                      "DPbla","CD4_SP","CD8_SP"))

#Set colors
colors <- c("DN3" = "#377EB8",
            "DN4" = "#4DAF4A",
            "DPel" = "#984EA3",
            "DPre" = "#FF7F00",
            "DPbla" = "#FFFF33",
            "CD4_SP" = "#A65628",
            "CD8_SP" = "#F781BF")

#Integrate data
DefaultAssay(scd.raw) <- 'RNA'
scd.raw[['SCT']] <- NULL
scd.raw[['integrated']] <- NULL

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
scd.cc <- scd.raw
scd.cc <- NormalizeData( scd.cc )
scd.cc <- FindVariableFeatures( scd.cc , selection.method = "vst" )
scd.cc <- ScaleData( scd.cc , features = rownames( scd.cc ) )
scd.cc <- CellCycleScoring(scd.cc, # min.cells = 1,
                           s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
cc.toadd <- data.frame( row.names = colnames(scd.cc) , S.Score = scd.cc$S.Score , G2M.Score = scd.cc$G2M.Score )
scd.raw <- AddMetaData( object = scd.raw , metadata = cc.toadd )
rm( scd.cc , cc.toadd )
cat('Cell Cycle Score complete.\n')

# Split Object and run SCTransformm while regressing out both mitochondrial genes and cell cycle genes
scd.list <- SplitObject(scd.raw, split.by = "Infection") # NORMALIZE BY SAMPLE #
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
scd <- IntegrateData(anchorset = scd.anchors, normalization.method = "SCT")
cat('Integration complete.\n')

scd <- RunPCA(scd, verbose = FALSE)
scd <- RunUMAP(scd, reduction = "pca", dims = 1:30)

scd <- FindNeighbors(scd, reduction = 'pca', dims = 1:30, verbose = FALSE)
scd <- FindClusters(scd)

##Save obj
saveRDS(scd, file.path(wd, "scd.sub.rds"))

gc()