#DEG analysis of DN
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
scd <- subset(scd.raw, (lin_X2 %in% c('DN1/2', 'DN3', 'DN4')))
table(scd$lin_X2) #Sanity check

Idents(scd) <- 'lin_X2'
DN <- scd
DN1.2. <- subset(scd, lin_X2 %in% 'DN1/2')
DN3 <- subset(scd, lin_X2 %in% 'DN3')
DN4<- subset(scd, lin_X2 %in% 'DN4')

DefaultAssay(DN) <- 'RNA'
DefaultAssay(DN1.2.) <- 'RNA'
DefaultAssay(DN3) <- 'RNA'
DefaultAssay(DN4) <- 'RNA'

DN <- ScaleData(DN)
DN1.2. <- ScaleData(DN1.2.)
DN3 <- ScaleData(DN3)
DN4<- ScaleData(DN4)

#DN
Idents(DN) <- 'Infection'
scd.markers <- FindAllMarkers(DN,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scd.markers <- subset(scd.markers, p_val_adj < 0.05)
scd.markers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> DN.markers
scd.Heatmap <- DoHeatmap(DN, features = DN.markers$gene, group.by = 'ident')

file.name <- paste('DN','Heatmap','pdf', sep='.')
pdf( file = file.path(wd, 'write', 'DN', file.name ) )
print( scd.Heatmap )
dev.off()

#DN1.2.
Idents(DN1.2.) <- 'Infection'
scd.markers <- FindAllMarkers(DN1.2.,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scd.markers <- subset(scd.markers, p_val_adj < 0.05)
scd.markers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> DN1.2..markers

#DN3
Idents(DN3) <- 'Infection'
scd.markers <- FindAllMarkers(DN3,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scd.markers <- subset(scd.markers, p_val_adj < 0.05)
scd.markers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> DN3.markers

#DN4
Idents(DN4) <- 'Infection'
scd.markers <- FindAllMarkers(DN4,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scd.markers <- subset(scd.markers, p_val_adj < 0.05)
scd.markers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> DN4.markers

##Pathway analysis
Mock.DN <- subset(DN.markers, cluster == 'Mock')
MRV.DN <- subset(DN.markers, cluster == 'MRV')
Mock.DN1.2. <- subset(DN1.2..markers, cluster == 'Mock')
MRV.DN1.2. <- subset(DN1.2..markers, cluster == 'MRV')
Mock.DN3 <- subset(DN3.markers, cluster == 'Mock')
MRV.DN3 <- subset(DN3.markers, cluster == 'MRV')
Mock.DN4 <- subset(DN4.markers, cluster == 'Mock')
MRV.DN4 <- subset(DN4.markers, cluster == 'MRV')

gene.list <- list(Mock.DN = Mock.DN,
                  MRV.DN = MRV.DN,
                  Mock.DN1.2. = Mock.DN1.2.,
                  MRV.DN1.2. = MRV.DN1.2.,
                  Mock.DN3 = Mock.DN3,
                  MRV.DN3 = MRV.DN3,
                  Mock.DN4 = Mock.DN4,
                  MRV.DN4 = MRV.DN4)

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
  plot <- if (websiteLive) plotEnrich(enriched[[1]], showTerms = g, numChar = 40, y = "Count", orderBy = "Adjusted.Adjusted.P.value",
                                      title = names(gene.list[i])) + 
    theme(text = element_text(size = 20))
  
  file.name <- paste( names(gene.list[i]), 'pdf', sep='.')
  pdf( file = file.path(wd, 'write', 'DN', 'EnrichR', file.name ) )
  print( plot )
  dev.off() 
  
  enriched[[1]] <- enriched[[1]] %>% 
    mutate(num = str_count(Genes, ";")+1)
  
  enriched[[1]] <- filter(enriched[[1]], Adjusted.P.value < 0.05)
  
  file.name <- paste( names(gene.list[i]), 'csv', sep='.')
  write.csv(enriched[[1]], file.path(wd, 'write', 'DN', 'Raw', file.name))
  
  names(enriched)[1] <- names(gene.list[i])
  csv.list[[i]] <- enriched[1]
}

####Combined EnrichR plots
#DN
DN.Mock <- data.frame(csv.list[1])
DN.MRV <- data.frame(csv.list[2])

colnames(DN.Mock) <- colnames(enriched[[1]])
colnames(DN.MRV) <- colnames(enriched[[1]])

DN.Mock$Term  <- gsub(" WP.*","", DN.Mock$Term)
DN.MRV$Term  <- gsub(" WP.*","", DN.MRV$Term)

DN.Mock$Term <- paste0(DN.Mock$Term, "-Mock")
DN.MRV$Term <- paste0(DN.MRV$Term, "-MRV")

#make MRV negative
DN.Mock$num <- -DN.Mock$num

#Reverse Mock df
DN.Mock%>% map_df(rev) -> DN.Mock

DN.both <- rbind(DN.MRV, DN.Mock)

#make plot
DN.both$Term <- factor(DN.both$Term, levels = rev(unique(DN.both$Term)))
plot <- ggplot(DN.both, aes(x = Term, y = num)) +
  geom_col(aes(fill = Adjusted.P.value), width = 0.7) + 
  coord_flip() +
  ggtitle("DN.Both") +
  ylab("Genes\nMock        MRV") +
  xlab("Pathway") +
  geom_hline(yintercept=0) +
  scale_fill_gradient(low="red",high="blue", trans = 'log', breaks=c(0.00005,0.0005,0.005,0.05),limits=c(0.000049,0.051)) +
  scale_y_continuous(breaks = round(seq(min(DN.both$num), max(DN.both$num), by = 4),1)) +
  theme_classic()

file.name <- paste( "DN.Both", 'pdf', sep='.')
pdf( file = file.path(wd, 'write', 'DN', 'Combined', file.name ) )
print( plot )
dev.off() 

#DN1
DN1.2.Mock <- data.frame(csv.list[3])
DN1.2.MRV <- data.frame(csv.list[4])

colnames(DN1.2.Mock) <- colnames(enriched[[1]])
colnames(DN1.2.MRV) <- colnames(enriched[[1]])

#DN1.2.Mock$Term  <- gsub(" WP.*","", DN1.2.Mock$Term)
DN1.2.MRV$Term  <- gsub(" WP.*","", DN1.2.MRV$Term)

#DN1.2.Mock$Term <- paste0(DN1.2.Mock$Term, "-Mock")
DN1.2.MRV$Term <- paste0(DN1.2.MRV$Term, "-MRV")

DN1.2.both <- DN1.2.MRV

#make plot
DN1.2.both$Term <- factor(DN1.2.both$Term, levels = rev(unique(DN1.2.both$Term)))
plot <- ggplot(DN1.2.both, aes(x = Term, y = num)) +
  geom_col(aes(fill = Adjusted.P.value), width = 0.7) + 
  coord_flip() +
  ggtitle("DN1.2.Both") +
  #ylab("Genes\nMock        MRV") +
  ylab("Genes\nMRV") +
  xlab("Pathway") +
  geom_hline(yintercept=0) +
  scale_fill_gradient(low="red",high="blue", trans = 'log', breaks=c(0.00005,0.0005,0.005,0.05),limits=c(0.000049,0.051)) +
  scale_y_continuous(breaks = round(seq(min(DN1.2.both$num), max(DN1.2.both$num), by = 4),1)) +
  theme_classic()

file.name <- paste( "DN1.2.Both", 'pdf', sep='.')
pdf( file = file.path(wd, 'write', 'DN', 'Combined', file.name ) )
print( plot )
dev.off() 

#DN3
DN3.Mock <- data.frame(csv.list[5])
DN3.MRV <- data.frame(csv.list[6])

colnames(DN3.Mock) <- colnames(enriched[[1]])
colnames(DN3.MRV) <- colnames(enriched[[1]])

DN3.Mock$Term  <- gsub(" WP.*","", DN3.Mock$Term)
DN3.MRV$Term  <- gsub(" WP.*","", DN3.MRV$Term)

DN3.Mock$Term <- paste0(DN3.Mock$Term, "-Mock")
DN3.MRV$Term <- paste0(DN3.MRV$Term, "-MRV")

#make MRV negative
DN3.Mock$num <- -DN3.Mock$num

#Reverse Mock df
DN3.Mock%>% map_df(rev) -> DN3.Mock

DN3.both <- rbind(DN3.MRV, DN3.Mock)

#make plot
DN3.both$Term <- factor(DN3.both$Term, levels = rev(unique(DN3.both$Term)))
plot <- ggplot(DN3.both, aes(x = Term, y = num)) +
  geom_col(aes(fill = Adjusted.P.value), width = 0.7) + 
  coord_flip() +
  ggtitle("DN3.Both") +
  ylab("Genes\nMock        MRV") +
  xlab("Pathway") +
  geom_hline(yintercept=0) +
  scale_fill_gradient(low="red",high="blue", trans = 'log', breaks=c(0.00005,0.0005,0.005,0.05),limits=c(0.000049,0.051)) +
  scale_y_continuous(breaks = round(seq(min(DN3.both$num), max(DN3.both$num), by = 4),1)) +
  theme_classic()

file.name <- paste( "DN3.Both", 'pdf', sep='.')
pdf( file = file.path(wd, 'write', 'DN', 'Combined', file.name ) )
print( plot )
dev.off() 

#DN4
DN4.Mock <- data.frame(csv.list[7])
DN4.MRV <- data.frame(csv.list[8])

colnames(DN4.Mock) <- colnames(enriched[[1]])
colnames(DN4.MRV) <- colnames(enriched[[1]])

#DN4.Mock$Term  <- gsub(" WP.*","", DN4.Mock$Term)
DN4.MRV$Term  <- gsub(" WP.*","", DN4.MRV$Term)

#DN4.Mock$Term <- paste0(DN4.Mock$Term, "-Mock")
DN4.MRV$Term <- paste0(DN4.MRV$Term, "-MRV")

DN4.both <- DN4.MRV

#make plot
DN4.both$Term <- factor(DN4.both$Term, levels = rev(unique(DN4.both$Term)))
plot <- ggplot(DN4.both, aes(x = Term, y = num)) +
  geom_col(aes(fill = Adjusted.P.value), width = 0.7) + 
  coord_flip() +
  ggtitle("DN4.Both") +
 # ylab("Genes\nMock        MRV") +
  ylab("Genes\nMRV") +
  xlab("Pathway") +
  geom_hline(yintercept=0) +
  scale_fill_gradient(low="red",high="blue", trans = 'log', breaks=c(0.00005,0.0005,0.005,0.05),limits=c(0.000049,0.051)) +
  scale_y_continuous(breaks = round(seq(min(DN4.both$num), max(DN4.both$num), by = 4),1)) +
  theme_classic()

file.name <- paste( "DN4.Both", 'pdf', sep='.')
pdf( file = file.path(wd, 'write', 'DN', 'Combined', file.name ) )
print( plot )
dev.off() 

#####Calculations
table <- table(scd$Infection, scd$lin_X2)
dft <- as.data.frame(rbind(table))
write.csv( table, file.path( wd, 'DN_Table.csv'))

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

file.name <- paste( 'DN_Bar','pdf', sep='.')
pdf( file = file.path(wd, 'write', 'DN', file.name ) )
print( plot )
dev.off()

#print calculations
write.csv(DN.markers, file.path(wd, 'write', 'DN', 'DN.csv'))
write.csv(DN1.2..markers, file.path(wd, 'write', 'DN', 'DN1.2..csv'))
write.csv(DN3.markers, file.path(wd, 'write', 'DN', 'DN3.csv'))
write.csv(DN4.markers, file.path(wd, 'write', 'DN', 'DN4.csv'))

gc()