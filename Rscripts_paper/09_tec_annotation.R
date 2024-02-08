#Annotate object and identify infected/replicating
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
options(future.globals.maxSize = 32000 * 1024^2) #32GB

# Load the Seurat Object 'scd.tec' ("Single Cell Data") or 'scd.tec' ("T cell subset")
wd <- getwd()
source( file = file.path(wd, 'src', 'shared_r_functions.R' ) )
rds.path <- file.path(wd, 'TEC.rds')
scd.tec <- readRDS(rds.path)
scd.tec # Sanity check
head( Idents(object = scd.tec) ) # Sanity check

# Read in Lineage Annotation File and Add Annotation as MetaData to Seurat Object than convert to Idents
Idents(scd.tec) <- "seurat_clusters"
lin.path <- file.path( annotation_path , "tec_annotation.tsv" )
lin_annotation <- read.delim( lin.path ) 
chosen_clusters <- data.frame( cluster = Idents(scd.tec) )
chosen_clusters$cluster <- as.numeric(as.character(chosen_clusters$cluster)) # This step was needed because the Idents were Factors that had Levels that did not correspond numerically
new.cluster.ids <- left_join( chosen_clusters, lin_annotation, by = "cluster" ) 
rownames(new.cluster.ids) <- rownames(chosen_clusters)
scd.tec <- AddMetaData( scd.tec, new.cluster.ids$lineage, col.name = 'TEC')

scd.tec$TEC <- factor(scd.tec$TEC, c("ETP","cTEC","mTEC",'DC',"Mono/Mac","Granulocyte","Endo/Mes"))

colors <- c("cTEC" = "red",
            'DC' = "steelblue",
            "Endo/Mes" = "green",
            "ETP" = "purple",
            "Granulocyte" = "orange",
            "Mono/Mac" = "yellow",
            "mTEC" = "brown")

Idents(scd.tec) <- "TEC"
plot <- UMAPPlot(scd.tec, split.by = "Infection", cols = colors)

file.name <- paste('TEC.Infection.annotation','pdf', sep='.')
pdf( file = file.path(wd, 'write.6', file.name ) )
print( plot )
dev.off()

#markers+Heatmap
scd.markers <- FindAllMarkers(scd.tec,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scd.markers <- subset(scd.markers, p_val_adj < 0.05)
scd.markers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> scd.tec.markers

write.csv(scd.tec.markers, file.path(wd, 'write.6', 'TEC.Markers.csv'))

#### Heatmap####
Idents(scd.tec) <- "Infection"
scd.markers <- FindAllMarkers(scd.tec,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scd.markers <- subset(scd.markers, p_val_adj < 0.05)
scd.markers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> scd.tec.markers

write.csv(scd.tec.markers, file.path(wd, 'write.6', 'TEC.Infection.Markers.csv'))

####define infection####
DefaultAssay(scd.tec) <- "virus"
inf <- subset(scd.tec, ORF55 >= 1 & ORF83 >= 1 & ORF103 >= 1)
rep <- subset(scd.tec,
              ( ORF53>=3 & ORF69>=3 & ORF73>=3 & ORF86>=3 ) |
                ( ORF45>=3 & ORF69>=3 & ORF73>=3 & ORF86>=3 ) |
                ( ORF45>=3 & ORF53>=3 & ORF73>=3 & ORF86>=3 ) |
                ( ORF45>=3 & ORF53>=3 & ORF69>=3 & ORF86>=3 ) |
                ( ORF45>=3 & ORF53>=3 & ORF69>=3 & ORF73>=3 ))

inf.barcodes <- data.frame(barcodes = rownames(inf@meta.data))
inf.barcodes$inf <- "Infected"
rep.barcodes <- data.frame(barcodes = rownames(rep@meta.data))
rep.barcodes$rep <- "Replicating"

#make lists for left_joins/coalesces
inf.list <- data.frame(barcodes = rownames(scd.tec@meta.data))
rep.list <- data.frame(barcodes = rownames(scd.tec@meta.data))
join.list <- data.frame(barcodes = rownames(scd.tec@meta.data))

#Identify overlap and differences between lists
dif.inf <- data.frame(barcodes = setdiff(inf.barcodes$barcodes, rep.barcodes$barcodes))
dif.rep <- data.frame(barcodes = setdiff(rep.barcodes$barcodes, inf.barcodes$barcodes))
same.both <- data.frame(barcodes = intersect(inf.barcodes$barcodes, rep.barcodes$barcodes))

print(length(dif.inf$barcodes) + length(same.both$barcodes))
print(length(inf.barcodes$barcodes))

print(length(dif.rep$barcodes) + length(same.both$barcodes))
print(length(rep.barcodes$barcodes))#Sanity checks, sets of values should be equal

#Add list data
dif.inf$Infected <- "Infected"
dif.rep$Infected <- "Replicating"
same.both$Infected <- "Replicating"

#left_join lists
inf.list <- left_join(join.list, dif.inf, "barcodes")
rep.list <- left_join(join.list, dif.rep, "barcodes")
join.list <- left_join(join.list, same.both, "barcodes")

print(table(inf.list$Infected))
print(table(rep.list$Infected))
print(table(join.list$Infected))#Sanity check

#Coalesce into one list
join.list$Infected <- coalesce(join.list$Infected,inf.list$Infected)
join.list$Infected <- coalesce(join.list$Infected,rep.list$Infected)

print(table(join.list$Infected))#Sanity check

#Add none values
join.list$Infected [is.na(join.list$Infected)] <- "None"

print(table(join.list$Infected))#Sanity check

##addmetadata
scd.tec <- AddMetaData(scd.tec, join.list$Infected, "Infected")
table(scd.tec$Infected)#Sanity check

n <- length(scd.tec$Infected)
for(i in 1:n){
  if(colSums(scd.tec@assays$virus@counts)[i] == 0){ #select value from plot
    scd.tec$Infected[i] <- "None"
  } else if(colSums(scd.tec@assays$virus@counts)[i] < 20){ # remove 0 values to make comparison
    scd.tec$Infected[i] <- "Un-Infected"
  }
}

###Heatmap
DefaultAssay(scd.tec) <- 'RNA'
Idents(scd.tec) <- "Infected"
scd.markers <- FindAllMarkers(scd.tec,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scd.markers <- subset(scd.markers, p_val_adj < 0.05)
scd.markers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> scd.tec.markers

write.csv(scd.tec.markers, file.path(wd, 'write.6', 'TEC.Infected.Markers.csv'))

#Infected percents
table <- as.data.frame.array(table(scd.tec$Infected, scd.tec$TEC))

write.csv(table, file.path(wd, 'write.6', 'TEC.Table.csv'))

percent <- sweep(table,2,colSums(table),`/`) * 100
percent.2 <- aggregate(percent, list(Group=replace(rownames(percent),rownames(percent) %in% c("Both","Replicating"), "Both+Replicating")), sum)
rownames(percent.2) <- percent.2$Group
percent.2$Group <- NULL

write.csv(percent.2, file.path(wd, 'write.6', 'TEC.Percent.csv'))


#Save file
saveRDS(scd.tec, file.path(wd, "TEC.rds"))

###Prop bar
table <- table(scd.tec$Infection, scd.tec$TEC)
dft <- as.data.frame(rbind(table))
write.csv( dft, file.path( wd, 'DP_Table.csv'))

#plot calculations
dft2 <- dft
dft2$Infection <- rownames(dft)
dft2 <- melt(dft2)
plot <- ggplot(dft2, aes(fill = variable, y = value, x = Infection)) +  scale_fill_manual(values = c(colors)) +  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous() +
  ylab("Proportion")

file.name <- paste( 'TEC.MRV.Mock_Bar','pdf', sep='.')
pdf( file = file.path(wd, 'write.6', file.name ) )
print( plot )
dev.off()


####ORF List####
orf.list <- data.frame(gene = c("ORF1",
                                "ORF2",
                                "ORF3",
                                "ORF4",
                                "ORF5",
                                "ORF6",
                                "ORF7",
                                "ORF8",
                                "ORF9",
                                "ORF10",
                                "ORF11",
                                "ORF12",
                                "ORF13",
                                "ORF14",
                                "ORF15",
                                "ORF16",
                                "ORF17",
                                "ORF18",
                                "ORF19",
                                "ORF20",
                                "ORF21",
                                "ORF22",
                                "ORF23",
                                "ORF24",
                                "ORF25",
                                "ORF26",
                                "ORF27",
                                "ORF28",
                                "ORF29",
                                "ORF30",
                                "ORF31",
                                "ORF32",
                                "ORF33",
                                "ORF34",
                                "ORF35",
                                "ORF36",
                                "ORF37",
                                "ORF38",
                                "ORF39",
                                "ORF40",
                                "ORF41",
                                "ORF42",
                                "ORF43",
                                "ORF44",
                                "ORF45",
                                "ORF46",
                                "ORF47",
                                "ORF48",
                                "ORF49",
                                "ORF50",
                                "ORF51",
                                "ORF52",
                                "ORF53",
                                "ORF54",
                                "ORF55",
                                "ORF56",
                                "ORF57",
                                "ORF58",
                                "ORF59",
                                "ORF60",
                                "ORF61",
                                "ORF62",
                                "ORF63",
                                "ORF64",
                                "ORF65",
                                "ORF66",
                                "ORF67",
                                "ORF68",
                                "ORF69",
                                "ORF70",
                                "ORF71",
                                "ORF72",
                                "ORF73",
                                "ORF74",
                                "ORF75",
                                "ORF76",
                                "ORF77",
                                "ORF78",
                                "ORF79",
                                "ORF80",
                                "ORF81",
                                "ORF82",
                                "ORF83",
                                "ORF84",
                                "ORF85",
                                "ORF86",
                                "ORF87",
                                "ORF88",
                                "ORF89",
                                "ORF90",
                                "ORF91",
                                "ORF92",
                                "ORF93",
                                "ORF94",
                                "ORF93",
                                "ORF96",
                                "ORF97",
                                "ORF98",
                                "ORF99",
                                "ORF100",
                                "ORF101",
                                "ORF102",
                                "ORF103",
                                "ORF104",
                                "ORF105",
                                "ORF106",
                                "ORF107",
                                "ORF108",
                                "ORF109",
                                "ORF110",
                                "ORF111",
                                "ORF112",
                                "ORF113",
                                "ORF114",
                                "ORF115",
                                "ORF116",
                                "ORF117",
                                "ORF118",
                                "ORF119",
                                "ORF120",
                                "ORF121",
                                "ORF122",
                                "ORF123",
                                "ORF124",
                                "ORF125",
                                "ORF126",
                                "ORF127",
                                "ORF128"))

orf.list$kinetics <- c("U",
                       "U",
                       "U",
                       "U",
                       "U",
                       "U",
                       "E/L",
                       "E/L",
                       "U",
                       "U",
                       "U",
                       "U",
                       "U",
                       "E/L",
                       "E/L",
                       "E/L",
                       "E/L",
                       "E/L",
                       "E/L",
                       "L",
                       "E",
                       "L",
                       "U",
                       "L",
                       "L",
                       "E",
                       "E/L",
                       "L",
                       "E",
                       "IE",
                       "IE",
                       "IE",
                       "L",
                       "L",
                       "U",
                       "U",
                       "L",
                       "U",
                       "E/L",
                       "E/L",
                       "E",
                       "L",
                       "U",
                       "U",
                       "E/L",
                       "U",
                       "E",
                       "L",
                       "L",
                       "E",
                       "E",
                       "U",
                       "E",
                       "E",
                       "E",
                       "E/L",
                       "E",
                       "E/L",
                       "E/L",
                       "U",
                       "L",
                       "E/L",
                       "L",
                       "E",
                       "E",
                       "E",
                       "L",
                       "U",
                       "L",
                       "U",
                       "E/L",
                       "E",
                       "E/L",
                       "E",
                       "U",
                       "U",
                       "L",
                       "L",
                       "L",
                       "E/L",
                       "U",
                       "U",
                       "E/L",
                       "U",
                       "E/L",
                       "L",
                       "U",
                       "E",
                       "L",
                       "U",
                       "E",
                       "L",
                       "U",
                       "U",
                       "E/L",
                       "U",
                       "U",
                       "E",
                       "L",
                       "U",
                       "U",
                       "E",
                       "IE",
                       "U",
                       "U",
                       "U",
                       "U",
                       "U",
                       "U",
                       "E",
                       "U",
                       "U",
                       "U",
                       "U",
                       "U",
                       "U",
                       "U",
                       "U",
                       "U",
                       "U",
                       "U",
                       "U",
                       "U",
                       "U",
                       "U",
                       "U",
                       "U",
                       "U") 
####ORF Heatmap ####
scd.tec <- subset(scd.tec, Infection == "MRV")
DefaultAssay(scd.tec) <- "virus"
orf.sub <- subset(orf.list, kinetics != "U")
table(orf.list$kinetics)
table(orf.sub$kinetics) #Sanity check

#reorder list by kinetics
scd.markers <- FindAllMarkers(scd.tec)

scd.markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) -> scd.markers

scd.markers <- subset(scd.markers, gene %in% orf.list$gene)
scd.markers <- scd.markers[match(orf.list$gene, scd.markers$gene),]

x <- c("IE","E","E/L","L")
orf.sub %>%
  mutate(kinetics =  factor(kinetics, levels = x)) %>%
  arrange(kinetics) -> orf.sub

scd.sub.markers <- subset(scd.markers, gene %in% orf.sub$gene)
scd.sub.markers <- scd.sub.markers[match(orf.sub$gene, scd.sub.markers$gene),]

scd.tec <- ScaleData(scd.tec)

##Multibar
scd.tec <- subset(scd.tec, Infected %in% c("Infected","Replicating","Un-Infected"))
scd.tec$Infected <- factor(scd.tec$Infected, c("Infected","Replicating","Un-Infected"))

cols.inf <- c("Un-Infected"="grey",
              "Infected"="red",
              "Replicating"="blue")
colors.to.use <- list( TEC = colors,  Infected = cols.inf)
scd.tec.Heatmap <- DoMultiBarHeatmap(scd.tec,
                                     features = scd.sub.markers$gene,
                                     group.by = 'TEC', additional.group.by = c('Infected'), cols.use = colors.to.use,
                                     additional.group.sort.by = 'Infected',
                                     angle = 90,
                                     size = 2,
                                     group.bar.height = 0.01,
                                     raster = F)  +
  theme(axis.text.y = element_text(size = 5), plot.margin = margin(40, 65, 0, 0)) + # decrease size of gene names
  # scale_fill_gradientn(colours = rev(mapal)) + 
  NoLegend()

file.name <- paste('ORF.sub.multi', 'pdf', sep='.')
pdf( file = file.path(wd, 'write.6', file.name ))
print( scd.tec.Heatmap )
dev.off()

gc()