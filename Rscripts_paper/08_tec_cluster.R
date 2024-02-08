#Make new dimensionality reduction for TECs/DCs
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
#library( DoMultiBarHeatmap )
# library( grid )
# library( escape )
library( VGAM )
library( patchwork )
library( ggpubr )
library( harmony )

#Increase memory usage
options(future.globals.maxSize = 32000 * 1024^2) #32GB

# Load the Seurat Object 'tsub' ("Single Cell Data") or 'tsub' ("T cell subset")
wd <- getwd()
rds.path <- file.path(wd, 'scd.lin_3.integrated.new.rds')
tsub <- readRDS(rds.path)
tsub # Sanity check
head( Idents(object = tsub) ) # Sanity check

#Remove assays
tsub@assays$SCT <- NULL
DefaultAssay(tsub) <- "RNA"

scd123 <- subset(tsub, lin_X2 %in% c("DC/ETP/Mac","Endo/Mes/cTEC/mTEC"))
scd123@assays$integrated<- NULL
DefaultAssay(scd123) <- "RNA"
####Process data for Harmony
###remove some variable sequences
#VariableFeatures(scd.list) <- VariableFeatures(scd.list)[!grepl("^TRAV|^TRAJ|^TRAC|^TRBV|^TRBD|^TRBJ|^TRBC|^IGHV|^IGHD|^IGHJ|^IGHM|^IGHE|^IGHG|^IGHA|^IGLV|^IGLD|^IGLC|^IGKV|^IGKD|^IGKC", scd.list@assays$RNA@var.features)]
counts <- GetAssayData(scd123, assay = "RNA")
vdj.genes <- rownames(counts)[grepl("^TRAV|^TRAJ|^TRBV|^TRBD|^TRBJ|^IGHV|^IGHD|^IGHJ|^IGLV|^IGLD|^IGKV|^IGKD", rownames(counts))]
counts <- counts[-(which(rownames(counts) %in% vdj.genes)),]
scd.p <- subset(scd123, features = rownames(counts))
scd.p <- NormalizeData(scd.p)
scd.p <- FindVariableFeatures(scd.p)
VariableFeatures(scd.p)[1:20]#sanity check
scd123 <- NormalizeData(scd123)
all.genes <- rownames(scd123)
scd123<- ScaleData(scd123, features = all.genes)
#Run PCA on filtered data
scd123 <- RunPCA(scd123, features = VariableFeatures(scd.p))
#Run Harmony on processed data


scd <- RunHarmony(scd123, group.by.vars = "Sample")
DefaultAssay(scd) <- "RNA"
#Run Umap
scd <- RunUMAP(scd, reduction = "harmony", dims = 1:30, verbose = FALSE) 

scd <- FindNeighbors(scd, reduction = 'harmony', dims = 1:30, verbose = FALSE, assay = "RNA")


scd <- FindClusters(scd, resolution = 2.0, verbose = FALSE)

#identify cluster based on markers, using umaps and violins to identify new clusters
Idents(scd) <- 'RNA_snn_res.2'
scd.markers <- FindAllMarkers(scd, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.csv(scd.markers, file.path(wd, "write.1", "All.mTEC.markers.csv"))

scd.markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < .05) %>%
  filter(pct.2 > .3) %>%
  filter(avg_log2FC > 1) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) -> scd.markers

saveRDS(scd, file.path(wd, "TEC.rds"))