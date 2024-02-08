#SP DEG analysis
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
options(future.globals.maxSize = 120000 * 1024^2) #120GB

# Load the Seurat Object 'tsub' ("Single Cell Data") or 'tsub' ("T cell subset")
wd <- getwd()
rds.path <- file.path(wd, 'scd.lin_3.integrated.new.rds')
scd.raw <- readRDS(rds.path)
scd.raw # Sanity check
head( Idents(object = scd.raw) ) # Sanity check

DefaultAssay(scd.raw) <- "RNA"

####CD vs DC/ETP/Mac(+ NK?) Heatmap
CD.sub <- subset(scd.raw, lin_X2 %in% c('CD4_SP', 'CD8_SP'))
#CD.sub <- subset(CD.sub, Cd4 == 0)
table(CD.sub$lin_X2)
table(scd.raw$lin_X2)#Sanity check


DefaultAssay(CD.sub) <- 'RNA'
CD.sub[['SCT']] <- NULL
CD.sub[['integrated']] <- NULL

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
scd.cc <- CD.sub
scd.cc <- NormalizeData( scd.cc )
scd.cc <- FindVariableFeatures( scd.cc , selection.method = "vst" )
scd.cc <- ScaleData( scd.cc , features = rownames( scd.cc ) )
scd.cc <- CellCycleScoring(scd.cc, # min.cells = 1,
                           s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
cc.toadd <- data.frame( row.names = colnames(scd.cc) , S.Score = scd.cc$S.Score , G2M.Score = scd.cc$G2M.Score )
CD.sub <- AddMetaData( object = CD.sub , metadata = cc.toadd )
rm( scd.cc , cc.toadd )
cat('Cell Cycle Score complete.\n')

# Split Object and run SCTransformm while regressing out both mitochondrial genes and cell cycle genes
scd.list <- SplitObject(CD.sub, split.by = "Infection") # NORMALIZE BY SAMPLE #
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
CD <- IntegrateData(anchorset = scd.anchors, normalization.method = "SCT", k.weight = 95)
cat('Integration complete.\n')

CD <- RunPCA(CD, verbose = FALSE)
CD <- RunUMAP(CD, reduction = "pca", dims = 1:30)

CD <- FindNeighbors(CD, reduction = 'pca', dims = 1:30, verbose = FALSE)
CD <- FindClusters(CD, resolution = .5)

# Print UMAP with chosen cluster resolution (located in meta.data under "integrated_snn_res.0.4" through "integrated_snn_res.0.6")
chosen_RES <- 'seurat_clusters'
#chosen_RES <- paste('integrated_snn_res', RES, sep = '.')
plot_clusters <- DimPlot(CD, group.by = chosen_RES, label = TRUE, pt.size = 0.025) + NoLegend()
file.name <- paste('Chosen' , 'clusters', '0.5', 'pdf', sep = '.')
pdf( file = file.path(wd, 'write.6', 'CD', file.name ) )
print( plot_clusters )
dev.off()

rm( file.name )
cat("Printed clustering.\n")

DefaultAssay(CD) <- "RNA"

####CD8/CD4 Degs
CD8 <- subset(CD, lin_X2 == 'CD8_SP')
CD4 <- subset(CD, lin_X2 == 'CD4_SP')

##CD8
Idents(CD8) <- 'Infection'
CDarkers <- FindAllMarkers(CD8,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CDarkers <- subset(CDarkers, p_val_adj < 0.05)
CDarkers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> CD8.markers
scd.Heatmap <- DoHeatmap(CD8, features = CD8.markers$gene, group.by = 'ident')

file.name <- paste('CD8','Heatmap','pdf', sep='.')
pdf( file = file.path(wd, 'write.6', 'CD8', file.name ) )
print( scd.Heatmap )
dev.off()

write.csv(CD8.markers, file.path(wd, 'write.6', 'CD8', 'CD8.markers.csv'))

##CD4
Idents(CD4) <- 'Infection'
CDarkers <- FindAllMarkers(CD4,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CDarkers <- subset(CDarkers, p_val_adj < 0.05)
CDarkers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> CD4.markers
scd.Heatmap <- DoHeatmap(CD4, features = CD4.markers$gene, group.by = 'ident')

file.name <- paste('CD4','Heatmap','pdf', sep='.')
pdf( file = file.path(wd, 'write.6', 'CD4', file.name ) )
print( scd.Heatmap )
dev.off()

write.csv(CD4.markers, file.path(wd, 'write.6', 'CD4', 'CD4.markers.csv'))

##EnrichR
Mock.CD8 <- subset(CD8.markers, cluster == 'Mock')
MRV.CD8 <- subset(CD8.markers, cluster == 'MRV')
Mock.CD4 <- subset(CD4.markers, cluster == 'Mock')
MRV.CD4 <- subset(CD4.markers, cluster == 'MRV')

gene.list <- list(Mock.CD8 = Mock.CD8,
                  MRV.CD8 = MRV.CD8,
                  Mock.CD4 = Mock.CD4,
                  MRV.CD4 = MRV.CD4)

n <- length(gene.list)

###General Enrichr pipeline
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes

websiteLive <- TRUE

dbs <- listEnrichrDbs()

if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

dbs <- c("WikiPathways_2019_Mouse", "HDSigDB_Mouse_2021", "KEGG_2019_Mouse")

#Function for wait
wait <- function(x)
{
  p1 <- proc.time()
  Sys.sleep(x)
  proc.time() - p1 # The cpu usage should be negligible
}

###Loop for to make pathway plots
csv.list <- list()
for(i in 1:n){
  wait(10)
  enriched <- c()
  if (websiteLive) {
    enriched <- enrichr(c(gene.list[[i]]$gene), dbs) #copy genes from heatmaps here
  }
  
  if (websiteLive) enriched[["WikiPathways_2019_Mouse"]]
  g <- nrow(filter(enriched[[1]], Adjusted.P.value < 0.05))
  plot <- if (websiteLive) plotEnrich(enriched[[1]], showTerms = g, numChar = 40, y = "Count", orderBy = "Adjusted.P.value",
                                      title = names(gene.list[i])) + 
    theme(text = element_text(size = 10))
  
  file.name <- paste( names(gene.list[i]), 'pdf', sep='.')
  pdf( file = file.path(wd, 'write.6', 'CD', 'EnrichR', file.name ) )
  print( plot )
  dev.off() 
  
  enriched[[1]] <- enriched[[1]] %>% 
    mutate(num = str_count(Genes, ";")+1)
  
  enriched[[1]] <- filter(enriched[[1]], Adjusted.P.value < 0.05)
  
  file.name <- paste( names(gene.list[i]), 'csv', sep='.')
  write.csv(enriched[[1]], file.path(wd, 'write.6', 'CD', 'Raw', file.name))
  
  names(enriched)[1] <- names(gene.list[i])
  csv.list[[i]] <- enriched[1]

}

####Combined EnrichR plots
#CD8
CD8.Mock <- data.frame(csv.list[1])
CD8.MRV <- data.frame(csv.list[2])

colnames(CD8.Mock) <- colnames(enriched[[1]])
colnames(CD8.MRV) <- colnames(enriched[[1]])

CD8.Mock$Term  <- gsub(" WP.*","", CD8.Mock$Term)
CD8.MRV$Term  <- gsub(" WP.*","", CD8.MRV$Term)

CD8.Mock$Term <- paste0(CD8.Mock$Term, "-Mock")
CD8.MRV$Term <- paste0(CD8.MRV$Term, "-MRV")

#make MRV negative
CD8.Mock$num <- -CD8.Mock$num

#Reverse Mock df
CD8.Mock%>% map_df(rev) -> CD8.Mock

CD8.both <- rbind(CD8.MRV, CD8.Mock)

#make plot
CD8.both$Term <- factor(CD8.both$Term, levels = rev(unique(CD8.both$Term)))
plot <- ggplot(CD8.both, aes(x = Term, y = num)) +
  geom_col(aes(fill = Adjusted.P.value), width = 0.7) + 
  coord_flip() +
  ggtitle("CD8.Both") +
  ylab("Genes\nMock        MRV") +
  xlab("Pathway") +
  geom_hline(yintercept=0) +
  scale_fill_gradient(low="red",high="blue", trans = 'log', breaks=c(0.00005,0.0005,0.005,0.05),limits=c(0.000049,0.051)) +
  scale_y_continuous(breaks = round(seq(min(CD8.both$num), max(CD8.both$num), by = 4),1)) +
  theme_classic()

file.name <- paste( "CD8.Both", 'pdf', sep='.')
pdf( file = file.path(wd, 'write.6', 'CD', 'Combined', file.name ) )
print( plot )
dev.off() 

#CD4
CD4.Mock <- data.frame(csv.list[3])
CD4.MRV <- data.frame(csv.list[4])

colnames(CD4.Mock) <- colnames(enriched[[1]])
colnames(CD4.MRV) <- colnames(enriched[[1]])

CD4.Mock$Term  <- gsub(" WP.*","", CD4.Mock$Term)
CD4.MRV$Term  <- gsub(" WP.*","", CD4.MRV$Term)

CD4.Mock$Term <- paste0(CD4.Mock$Term, "-Mock")
CD4.MRV$Term <- paste0(CD4.MRV$Term, "-MRV")

#make MRV negative
CD4.Mock$num <- -CD4.Mock$num

#Reverse Mock df
CD4.Mock%>% map_df(rev) -> CD4.Mock

CD4.both <- rbind(CD4.MRV, CD4.Mock)

#make plot
CD4.both$Term <- factor(CD4.both$Term, levels = rev(unique(CD4.both$Term)))
plot <- ggplot(CD4.both, aes(x = Term, y = num)) +
  geom_col(aes(fill = Adjusted.P.value), width = 0.7) + 
  coord_flip() +
  ggtitle("CD4.Both") +
  ylab("Genes\nMock        MRV") +
  xlab("Pathway") +
  geom_hline(yintercept=0) +
  scale_fill_gradient(low="red",high="blue", trans = 'log', breaks=c(0.00005,0.0005,0.005,0.05),limits=c(0.000049,0.051)) +
  scale_y_continuous(breaks = round(seq(min(CD4.both$num), max(CD4.both$num), by = 4),1)) +
  theme_classic()

file.name <- paste( "CD4.Both", 'pdf', sep='.')
pdf( file = file.path(wd, 'write.6', 'CD', 'Combined', file.name ) )
print( plot )
dev.off()

#Save RDS
saveRDS(CD4, file.path(wd, "CD4.rds"))
saveRDS(CD8, file.path(wd, "CD8.rds"))

gc()