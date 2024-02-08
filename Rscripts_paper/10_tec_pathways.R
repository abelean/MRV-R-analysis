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
library( VGAM )
library( patchwork )
library( ggpubr )
library( enrichR )
library(reshape2) 

#Increase memory usage
options(future.globals.maxSize = 120000 * 1024^2) #120GB

# Load the Seurat Object 'scd.tec' ("Single Cell Data") or 'scd.tec' ("T cell subset")
wd <- getwd()
source( file = file.path(wd, 'src', 'shared_r_functions.R' ) )
rds.path <- file.path(wd, 'TEC.rds')
scd.tec <- readRDS(rds.path)
scd.tec # Sanity check
head( Idents(object = scd.tec) ) # Sanity check

#Subset clusters
Mono.Mac <- subset(scd.tec, TEC == "Mono/Mac")
DC <- subset(scd.tec, TEC == "DC")
ETP <- subset(scd.tec, TEC == "ETP")
cTEC <- subset(scd.tec, TEC == "cTEC")
Granulocyte <- subset(scd.tec, TEC == "Granulocyte")
mTEC <- subset(scd.tec, TEC == "mTEC")
Endo.Mes <- subset(scd.tec, TEC == "Endo/Mes")

###Findallmarkers
object.list <- list(Mono.Mac = Mono.Mac,
                  DC = DC,
                  ETP = ETP,
                  cTEC = cTEC,
                  Granulocyte = Granulocyte,
                  mTEC = mTEC,
                  Endo.Mes = Endo.Mes)

scd.tec.markers <- list("Mono.Mac",
                        "DC",
                        "ETP",
                        "cTEC",
                        "Granulocyte",
                        "mTEC",
                        "Endo.Mes")

names(scd.tec.markers) <- list("Mono.Mac",
                               "DC",
                               "ETP",
                               "cTEC",
                               "Granulocyte",
                               "mTEC",
                               "Endo.Mes")

for(i in 1:7){
Idents(object.list[[i]]) <- "Infection"
scd.markers <- FindAllMarkers(object.list[[i]],  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scd.markers <- subset(scd.markers, p_val_adj < 0.05)
scd.markers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> scd.tec.markers[[i]]

if(length(scd.tec.markers[[i]]$p_val) > 0){
plot <- DoHeatmap(object.list[[i]], scd.tec.markers[[i]]$gene, group.by = 'ident', raster = F) +
  theme(axis.text.y = element_text(size = 8))

file.name <- paste('Heatmap', names(scd.tec.markers[i]), 'pdf', sep='.')
pdf( file = file.path(wd, 'write.6', 'EnrichR', file.name ) )
print( plot )
dev.off()

file.name <- paste('(Mock vs MRV)Markers', names(scd.tec.markers[i]), 'csv', sep='.')
write.csv(scd.tec.markers[[i]], file.path(wd,'write.6','EnrichR',file.name))
}
}

##Pathway analysis
Mock.Mono.Mac <- subset(scd.tec.markers[[1]], cluster == "Mock")
Mock.DC <- subset(scd.tec.markers[[2]], cluster == "Mock")
Mock.ETP <- subset(scd.tec.markers[[3]], cluster == "Mock")
Mock.cTEC <- subset(scd.tec.markers[[4]], cluster == "Mock")
Mock.Granulocyte <- subset(scd.tec.markers[[5]], cluster == "Mock")
Mock.mTEC <- subset(scd.tec.markers[[6]], cluster == "Mock")
Mock.Endo.Mes <- subset(scd.tec.markers[[7]], cluster == "Mock")

MRV.Mono.Mac <- subset(scd.tec.markers[[1]], cluster == "MRV")
MRV.DC <- subset(scd.tec.markers[[2]], cluster == "MRV")
MRV.ETP <- subset(scd.tec.markers[[3]], cluster == "MRV")
MRV.cTEC <- subset(scd.tec.markers[[4]], cluster == "MRV")
MRV.Granulocyte <- subset(scd.tec.markers[[5]], cluster == "MRV")
MRV.mTEC <- subset(scd.tec.markers[[6]], cluster == "MRV")
MRV.Endo.Mes <- subset(scd.tec.markers[[7]], cluster == "MRV")

gene.list <- list(Mock.Mono.Mac = Mock.Mono.Mac,
                  Mock.DC = Mock.DC,
                  Mock.ETP = Mock.ETP,
                  Mock.cTEC = Mock.cTEC,
                  Mock.Granulocyte = Mock.Granulocyte,
                  Mock.mTEC = Mock.mTEC,
                  Mock.Endo.Mes = Mock.Endo.Mes,
                  
                  MRV.Mono.Mac = MRV.Mono.Mac,
                  MRV.DC = MRV.DC,
                  MRV.ETP = MRV.ETP,
                  MRV.cTEC = MRV.cTEC,
                  MRV.Granulocyte = MRV.Granulocyte,
                  MRV.mTEC = MRV.mTEC,
                  MRV.Endo.Mes = MRV.Endo.Mes)

#overall markers
Idents(scd.tec) <- "TEC"
scd.markers <- FindAllMarkers(scd.tec,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scd.markers <- subset(scd.markers, p_val_adj < 0.05)
scd.markers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> scd.tec.markers.all

write.csv(scd.tec.markers.all, file.path(wd, 'write.6', 'TEC.markers.csv'))


scd.tec.markers <- list("Mono.Mac",
                        "DC",
                        "ETP",
                        "cTEC",
                        "Granulocyte",
                        "mTEC",
                        "Endo.Mes")

names(scd.tec.markers) <- list("Mono.Mac",
                               "DC",
                               "ETP",
                               "cTEC",
                               "Granulocyte",
                               "mTEC",
                               "Endo.Mes")


scd.tec.markers[[1]] <- subset(scd.tec.markers.all, cluster == "Mono/Mac")
scd.tec.markers[[2]] <- subset(scd.tec.markers.all, cluster == "DC")
scd.tec.markers[[3]] <- subset(scd.tec.markers.all, cluster == "ETP")
scd.tec.markers[[4]] <- subset(scd.tec.markers.all, cluster == "cTEC")
scd.tec.markers[[5]] <- subset(scd.tec.markers.all, cluster == "Granulocyte")
scd.tec.markers[[6]] <- subset(scd.tec.markers.all, cluster == "mTEC")
scd.tec.markers[[7]] <- subset(scd.tec.markers.all, cluster == "Endo/Mes")

x <- length(scd.tec.markers)

#Remove 0 length lists
for(i in 1:length(gene.list)){
  if(length(gene.list[[i]]$gene) == 0){
    gene.list[i] <- NA
  }
}

gene.list <- gene.list[!is.na(gene.list)]

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
for(i in 1:x){
  wait(10)
  enriched <- c()
  if (websiteLive) {
    enriched <- enrichr(c(scd.tec.markers[[i]]$gene), dbs) #copy genes from heatmaps here
  }
  
  if (websiteLive) enriched[["WikiPathways_2019_Mouse"]]
  
  enriched[[1]] <- enriched[[1]] %>% 
    mutate(num = str_count(Genes, ";")+1)
  
  enriched[[1]] <- filter(enriched[[1]], Adjusted.P.value < 0.05)
  enriched[[1]] %>% 
    arrange(Adjusted.P.value) -> enriched[[1]]
  
  file.name <- paste( names(scd.tec.markers[i]), 'csv', sep='.')
  write.csv(enriched[[1]], file.path(wd, 'write.6', 'Raw', file.name))

}

csv.list <- list()
for(i in 1:n){
  wait(10)
  enriched <- c()
  if (websiteLive) {
    enriched <- enrichr(c(gene.list[[i]]$gene), dbs) #copy genes from heatmaps here
  }
  
  if (websiteLive) enriched[["WikiPathways_2019_Mouse"]]
  
  enriched[[1]] <- enriched[[1]] %>% 
    mutate(num = str_count(Genes, ";")+1)
  
  enriched[[1]] <- filter(enriched[[1]], Adjusted.P.value < 0.05)
  enriched[[1]] %>% 
    arrange(Adjusted.P.value) -> enriched[[1]]
  
  file.name <- paste( names(gene.list[i]), 'csv', sep='.')
  write.csv(enriched[[1]], file.path(wd, 'write.6', 'Raw', file.name))

}

gc()