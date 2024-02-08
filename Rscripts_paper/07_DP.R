#DEG analysis of DP
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
library(scales) ###

#Increase memory usage
options(future.globals.maxSize = 32000 * 1024^2) #32GB

# Load the Seurat Object 'tsub' ("Single Cell Data") or 'tsub' ("T cell subset")
wd <- getwd()
rds.path <- file.path(wd, 'scd.lin_3.integrated.new.rds')
scd.raw <- readRDS(rds.path)
scd.raw # Sanity check
head( Idents(object = scd.raw) ) # Sanity check

####make heatmaps
Idents(scd.raw) <- 'lin_X2'
scd <- subset(scd.raw, (lin_X2 %in% c('DPbla', 'DPre', 'DPsel')))
table(scd$lin_X2) #Sanity check

Idents(scd) <- 'lin_X2'
DP <- scd
DPbla <- subset(scd, lin_X2 %in% 'DPbla')
DPre <- subset(scd, lin_X2 %in% 'DPre')
DPsel<- subset(scd, lin_X2 %in% 'DPsel')

DefaultAssay(DP) <- 'RNA'
DefaultAssay(DPbla) <- 'RNA'
DefaultAssay(DPre) <- 'RNA'
DefaultAssay(DPsel) <- 'RNA'

DP <- ScaleData(DP)
DPbla <- ScaleData(DPbla)
DPre <- ScaleData(DPre)
DPsel<- ScaleData(DPsel)

#DP
Idents(DP) <- 'Infection'
scd.markers <- FindAllMarkers(DP,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scd.markers <- subset(scd.markers, p_val_adj < 0.05)
scd.markers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> DP.markers
scd.Heatmap <- DoHeatmap(DP, features = DP.markers$gene, group.by = 'ident')

file.name <- paste('DP','Heatmap','pdf', sep='.')
pdf( file = file.path(wd, 'write', 'DP', file.name ) )
print( scd.Heatmap )
dev.off()

#DPbla
Idents(DPbla) <- 'Infection'
scd.markers <- FindAllMarkers(DPbla,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scd.markers <- subset(scd.markers, p_val_adj < 0.05)
scd.markers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> DPbla.markers

#DPre
Idents(DPre) <- 'Infection'
scd.markers <- FindAllMarkers(DPre,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scd.markers <- subset(scd.markers, p_val_adj < 0.05)
scd.markers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> DPre.markers

#DPsel
Idents(DPsel) <- 'Infection'
scd.markers <- FindAllMarkers(DPsel,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scd.markers <- subset(scd.markers, p_val_adj < 0.05)
scd.markers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> DPsel.markers

##Pathway analysis
Mock.DP <- subset(DP.markers, cluster == 'Mock')
MRV.DP <- subset(DP.markers, cluster == 'MRV')
Mock.DPbla <- subset(DPbla.markers, cluster == 'Mock')
MRV.DPbla <- subset(DPbla.markers, cluster == 'MRV')
Mock.DPre <- subset(DPre.markers, cluster == 'Mock')
MRV.DPre <- subset(DPre.markers, cluster == 'MRV')
Mock.DPsel <- subset(DPsel.markers, cluster == 'Mock')
MRV.DPsel <- subset(DPsel.markers, cluster == 'MRV')

gene.list <- list(Mock.DP = Mock.DP,
                  MRV.DP = MRV.DP,
                  Mock.DPbla = Mock.DPbla,
                  MRV.DPbla = MRV.DPbla,
                  Mock.DPre = Mock.DPre,
                  MRV.DPre = MRV.DPre,
                  Mock.DPsel = Mock.DPsel,
                  MRV.DPsel = MRV.DPsel)

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
    theme(text = element_text(size = 20))
  
  file.name <- paste( names(gene.list[i]), 'pdf', sep='.')
  pdf( file = file.path(wd, 'write', 'DP', 'EnrichR', file.name ) )
  print( plot )
  dev.off() 
  
  enriched[[1]] <- enriched[[1]] %>% 
    mutate(num = str_count(Genes, ";")+1)
  
  enriched[[1]] <- filter(enriched[[1]], Adjusted.P.value < 0.05)
  
  file.name <- paste( names(gene.list[i]), 'csv', sep='.')
  write.csv(enriched[[1]], file.path(wd, 'write', 'DP', 'Raw', file.name))
  
  names(enriched)[1] <- names(gene.list[i])
  csv.list[[i]] <- enriched[1]
}

####Combined EnrichR plots
#DP
DP.Mock <- data.frame(csv.list[1])
DP.MRV <- data.frame(csv.list[2])

colnames(DP.Mock) <- colnames(enriched[[1]])
colnames(DP.MRV) <- colnames(enriched[[1]])

DP.Mock$Term  <- gsub(" WP.*","", DP.Mock$Term)
DP.MRV$Term  <- gsub(" WP.*","", DP.MRV$Term)

DP.Mock$Term <- paste0(DP.Mock$Term, "-Mock")
DP.MRV$Term <- paste0(DP.MRV$Term, "-MRV")

#make MRV negative
DP.Mock$num <- -DP.Mock$num

#Reverse Mock df
DP.Mock%>% map_df(rev) -> DP.Mock

DP.both <- rbind(DP.MRV, DP.Mock)

#make plot
DP.both$Term <- factor(DP.both$Term, levels = rev(unique(DP.both$Term)))
plot <- ggplot(DP.both, aes(x = Term, y = num)) +
  geom_col(aes(fill = Adjusted.P.value), width = 0.7) + 
  coord_flip() +
  ggtitle("DP.Both") +
  ylab("Genes\nMock        MRV") +
  xlab("Pathway") +
  geom_hline(yintercept=0) +
  scale_fill_gradient(low="red",high="blue", trans = 'log') +
  scale_y_continuous(breaks = round(seq(min(DP.both$num), max(DP.both$num), by = 4),1)) +
  theme_classic()

file.name <- paste( "DP.Both", 'pdf', sep='.')
pdf( file = file.path(wd, 'write', 'DP', 'Combined', file.name ) )
print( plot )
dev.off() 

#DPbla
DPbla.Mock <- data.frame(csv.list[3])
DPbla.MRV <- data.frame(csv.list[4])

colnames(DPbla.Mock) <- colnames(enriched[[1]])
colnames(DPbla.MRV) <- colnames(enriched[[1]])

DPbla.Mock$Term  <- gsub(" WP.*","", DPbla.Mock$Term)
DPbla.MRV$Term  <- gsub(" WP.*","", DPbla.MRV$Term)

DPbla.Mock$Term <- paste0(DPbla.Mock$Term, "-Mock")
DPbla.MRV$Term <- paste0(DPbla.MRV$Term, "-MRV")

#make MRV negative
DPbla.Mock$num <- -DPbla.Mock$num

#Reverse Mock df
DPbla.Mock%>% map_df(rev) -> DPbla.Mock

DPbla.both <- rbind(DPbla.MRV, DPbla.Mock)

#make plot
DPbla.both$Term <- factor(DPbla.both$Term, levels = rev(unique(DPbla.both$Term)))
plot <- ggplot(DPbla.both, aes(x = Term, y = num)) +
  geom_col(aes(fill = Adjusted.P.value), width = 0.7) + 
  coord_flip() +
  ggtitle("DPbla.Both") +
  ylab("Genes\nMock        MRV") +
  xlab("Pathway") +
  geom_hline(yintercept=0) +
  scale_fill_gradient(low="red",high="blue", trans = 'log') +
  scale_y_continuous(breaks = round(seq(min(DPbla.both$num), max(DPbla.both$num), by = 4),1)) +
  theme_classic()

file.name <- paste( "DPbla.Both", 'pdf', sep='.')
pdf( file = file.path(wd, 'write', 'DP', 'Combined', file.name ) )
print( plot )
dev.off() 

#DPre
DPre.Mock <- data.frame(csv.list[5])
DPre.MRV <- data.frame(csv.list[6])

colnames(DPre.Mock) <- colnames(enriched[[1]])
colnames(DPre.MRV) <- colnames(enriched[[1]])

DPre.Mock$Term  <- gsub(" WP.*","", DPre.Mock$Term)
DPre.MRV$Term  <- gsub(" WP.*","", DPre.MRV$Term)

DPre.Mock$Term <- paste0(DPre.Mock$Term, "-Mock")
DPre.MRV$Term <- paste0(DPre.MRV$Term, "-MRV")

#make MRV negative
DPre.Mock$num <- -DPre.Mock$num

#Reverse Mock df
DPre.Mock%>% map_df(rev) -> DPre.Mock

DPre.both <- rbind(DPre.MRV, DPre.Mock)

#make plot
DPre.both$Term <- factor(DPre.both$Term, levels = rev(unique(DPre.both$Term)))
plot <- ggplot(DPre.both, aes(x = Term, y = num)) +
  geom_col(aes(fill = Adjusted.P.value), width = 0.7) + 
  coord_flip() +
  ggtitle("DPre.Both") +
  ylab("Genes\nMock        MRV") +
  xlab("Pathway") +
  geom_hline(yintercept=0) +
  scale_fill_gradient(low="red",high="blue", trans = 'log') +
  scale_y_continuous(breaks = round(seq(min(DPre.both$num), max(DPre.both$num), by = 4),1)) +
  theme_classic()

file.name <- paste( "DPre.Both", 'pdf', sep='.')
pdf( file = file.path(wd, 'write', 'DP', 'Combined', file.name ) )
print( plot )
dev.off() 

#DPsel
DPsel.Mock <- data.frame(csv.list[7])
DPsel.MRV <- data.frame(csv.list[8])

colnames(DPsel.Mock) <- colnames(enriched[[1]])
colnames(DPsel.MRV) <- colnames(enriched[[1]])

DPsel.Mock$Term  <- gsub(" WP.*","", DPsel.Mock$Term)
DPsel.MRV$Term  <- gsub(" WP.*","", DPsel.MRV$Term)

DPsel.Mock$Term <- paste0(DPsel.Mock$Term, "-Mock")
DPsel.MRV$Term <- paste0(DPsel.MRV$Term, "-MRV")

#make MRV negative
DPsel.Mock$num <- -DPsel.Mock$num

#Reverse Mock df
DPsel.Mock%>% map_df(rev) -> DPsel.Mock

DPsel.both <- rbind(DPsel.MRV, DPsel.Mock)

#make plot
DPsel.both$Term <- factor(DPsel.both$Term, levels = rev(unique(DPsel.both$Term)))
plot <- ggplot(DPsel.both, aes(x = Term, y = num)) +
  geom_col(aes(fill = Adjusted.P.value), width = 0.7) + 
  coord_flip() +
  ggtitle("DPsel.Both") +
  ylab("Genes\nMock        MRV") +
  xlab("Pathway") +
  geom_hline(yintercept=0) +
  scale_fill_gradient(low="red",high="blue", trans = 'log') +
  scale_y_continuous(breaks = round(seq(min(DPsel.both$num), max(DPsel.both$num), by = 4),1)) +
  theme_classic()

file.name <- paste( "DPsel.Both", 'pdf', sep='.')
pdf( file = file.path(wd, 'write', 'DP', 'Combined', file.name ) )
print( plot )
dev.off() 

#####Calculations
table <- table(scd$Infection, scd$lin_X2)
dft <- as.data.frame(rbind(table))
write.csv( table, file.path( wd, 'DP_Table.csv'))

#plot calculations
dft2 <- dft
dft2$Infection <- rownames(dft)
dft2 <- melt(dft2)
plot <- ggplot(dft2, aes(fill = variable, y = value, x = Infection)) + geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous() +
  scale_fill_manual(values = c("DN1/2" = "red",
                               'DN3' = "steelblue",
                               "DN4" = "green",
                               "DPbla" = "purple",
                               "DPre" = "orange",
                               "DPsel" = "yellow",
                               "CD4_SP" = "brown",
                               "CD8_SP" = "pink",
                               'ISP' = "turquoise",
                               'NK/ILC' = "magenta",
                               'Endo/Mes/cTEC/mTEC' = "beige",
                               'DC/ETP/Mac' = "violet",
                               'RBC' = "darkred")) +
  theme_bw() +
  ylab("Proportion")

file.name <- paste( 'DP_Bar','pdf', sep='.')
pdf( file = file.path(wd, 'write', 'DP', file.name ) )
print( plot )
dev.off()

#print calculations
write.csv(DP.markers, file.path(wd, 'write', 'DP', 'DP.csv'))
write.csv(DPbla.markers, file.path(wd, 'write', 'DP', 'DPbla.csv'))
write.csv(DPre.markers, file.path(wd, 'write', 'DP', 'DPre.csv'))
write.csv(DPsel.markers, file.path(wd, 'write', 'DP', 'DPsel.csv'))

gc()