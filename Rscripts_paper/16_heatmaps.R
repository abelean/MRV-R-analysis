#Generating heatmaps
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
#library( grid )
# library( escape )
#library( VGAM )
library( patchwork )
library( ggpubr )
library( harmony )
library( pheatmap )

#Increase memory usage
options(future.globals.maxSize = 32000 * 1024^2) #32GB

# Load the Seurat Object 'scd' ("Single Cell Data") or 'scd' ("T cell subset")
wd <- getwd()
source( file = file.path(wd, "src", "shared_r_functions.R" ) )

rds.path <- file.path(wd, 'scd.lin_3.integrated.new.rds')
scd <- readRDS(rds.path)
scd # Sanity check
head( Idents(object = scd) ) # Sanity check

# Load the Seurat Object 'scd' ("Single Cell Data")
rds.path <- file.path(wd, 'Infected.rds')
scd.inf <- readRDS(rds.path)
scd.inf # Sanity check
head( Idents(object = scd.inf) ) # Sanity check

#### Markers ####
Idents(scd) <- "Infection"
scd.markers <- FindAllMarkers(scd,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scd.markers <- subset(scd.markers, p_val_adj < 0.05)
scd.markers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> scd.markers

write.csv(scd.markers, file.path(wd,'write.2','(Mock vs MRV)All.markers.csv'))

Idents(scd) <- "seurat_clusters"
scd.markers <- FindAllMarkers(scd,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scd.markers <- subset(scd.markers, p_val_adj < 0.05)
scd.markers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> scd.markers

write.csv(scd.markers, file.path(wd,'write.2','(clusters)All.markers.csv'))

#subset out MRV
scd <- subset(scd, Infection == 'MRV')

####make list of orfs####
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
#### ####
#Apply labels to main object
barcodes <- data.frame("barcodes"=rownames(scd@meta.data))
infection <- data.frame("barcodes"=rownames(scd.inf@meta.data), "Infected" = scd.inf$Infected)
list.join <- left_join(barcodes, infection, by = "barcodes")

scd <- AddMetaData(scd, list.join$Infected, "Infected")

print(table(scd$Infected))
print(table(scd.inf$Infected))#Sanity check

#change none to uninfected
scd$Infected[is.na(scd$Infected)] <- 'Un-Infected'
scd$Infected[scd$Infected == 'None'] <- 'Un-Infected'
print(table(scd$Infected))#Sanity check

##make multibar heatmap
#markers
DefaultAssay(scd) <- "virus"
my_levels <- c('DN1/2', 'DN3', 'DN4', 'ISP', 'DPbla', 'DPre', 'DPsel', 'CD4_SP', 'CD8_SP', 'NK/ILC', 'Endo/Mes/cTEC/mTEC', 'DC/ETP/Mac', 'RBC')
scd$lin_X2 <- factor(scd$lin_X2, levels = c(my_levels))
scd$Infected <- factor(scd$Infected, levels =  c("Infected","Replicating","Un-Infected"))

Idents(scd) <- 'lin_X2'
scd <- ScaleData(scd)

scd.markers <- FindAllMarkers(scd)

scd.markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) -> scd.markers

scd.markers <- subset(scd.markers, gene %in% orf.list$gene)
scd.markers <- scd.markers[match(orf.list$gene, scd.markers$gene),]


orf.sub <- subset(orf.list, kinetics != "U")
table(orf.list$kinetics)
table(orf.sub$kinetics) #Sanity check

x <- c("IE","E","E/L","L")
orf.sub %>%
  mutate(kinetics =  factor(kinetics, levels = x)) %>%
  arrange(kinetics) -> orf.sub


scd.sub.markers <- subset(scd.markers, gene %in% orf.sub$gene)
scd.sub.markers <- scd.sub.markers[match(orf.sub$gene, scd.sub.markers$gene),]
#Heatmap
mapal <- colorRampPalette( RColorBrewer::brewer.pal( 11, "RdBu" ) )( 256 )
colors <- c("DN1/2" = "red",
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
            'RBC' = "darkred")
cols.inf <- c("Un-Infected"="grey",
              "Infected"="red",
              "Replicating"="blue")
colors.to.use <- list( lin_X2 = colors,  Infected = cols.inf)
scd.Heatmap <- DoMultiBarHeatmap(scd,
                                 features = scd.sub.markers$gene,
                                 group.by = 'lin_X2', additional.group.by = c('Infected'), cols.use = colors.to.use,
                                 additional.group.sort.by = 'Infected',
                                 angle = 90,
                                 size = 2,
                                 group.bar.height = 0.01,
                                 raster = F)  +
  theme(axis.text.y = element_text(size = 5), plot.margin = margin(40, 65, 0, 0)) + # decrease size of gene names
  # scale_fill_gradientn(colours = rev(mapal)) + 
  NoLegend()

file.name <- paste('ORF.sub.multi', 'pdf', sep='.')
pdf( file = file.path(wd, 'write.2', file.name ))
print( scd.Heatmap )
dev.off()

##ORF ordered heatmap 
Idents(scd) <- "Infected"

scd.markers <- FindAllMarkers(scd)
scd.markers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> scd.markers

Idents(scd) <- "lin_X2"

mapal <- colorRampPalette( RColorBrewer::brewer.pal( 11, "RdBu" ) )( 256 )
cols.inf <- c("Un-Infected"="grey",
              "Infected"="red",
              "Replicating"="blue")
colors.to.use <- list( lin_X2 = colors,  Infected = cols.inf)

#### Markers(cont.) ####
DefaultAssay(scd) <- "RNA"
Idents(scd) <- "Infected"
scd.markers <- FindAllMarkers(scd,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scd.markers <- subset(scd.markers, p_val_adj < 0.05)
scd.markers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> scd.markers

write.csv(scd.markers, file.path(wd,'write.2','(Infected)All.markers.csv'))

DefaultAssay(scd) <- "virus"

##Pheatmap
pheat <- data.frame(scd@assays$virus@counts)
pheat <-  data.frame(NormalizeData(pheat))

#Add other U
join.list <- data.frame("gene" = rownames(pheat))
join.list <- left_join(join.list, orf.list, "gene")
join.list$kinetics[is.na(join.list$kinetics)] <- "U"

#order rows
x <- c("IE","E","E/L","L","U")
join.list %>%
  mutate(kinetics =  factor(kinetics, levels = x)) %>%
  arrange(kinetics) -> join.list

pheat <- pheat[match(join.list$gene, rownames(pheat)), ]

print(pheat[1:10,1:10])
pheat <- pheat[, match(rownames(scd@meta.data), colnames(pheat))]
print(pheat[1:10,1:10])##Sanity check

sub_anno <- structure(list(seq_share = c(join.list$kinetics)), .Names = "Kinetics", row.names = c(rownames(pheat) ), class = "data.frame")

sub_anno2 <- structure(list(seq_share = c(scd$Infected)), .Names = "Infected", row.names = c(colnames(pheat) ), class = "data.frame")

ann_colors = list(Infected = c("Un-Infected"="grey",
                               "Infected"="red",
                               "Replicating"="blue"))

plot <- pheatmap(pheat,
         show_colnames = F,
         annotation_row = sub_anno,
         annotation_col = sub_anno2,
         fontsize_row = 5,
         annotation_colors = ann_colors)

file.name <- paste('ORF.p.heatmap', 'pdf', sep='.')
pdf( file = file.path(wd, 'write.2', file.name ), height = 10, width = 5)
print( plot )
dev.off()

####RNA Heatmaps
DefaultAssay(scd) <- "RNA"
scd <- ScaleData(scd)

scd.DPre <- subset(scd, lin_X2 == "DPre")
scd.DPsel <- subset(scd, lin_X2 == "DPsel")
scd.CD4SP <- subset(scd, lin_X2 == "CD4_SP")

#save
saveRDS(scd, file.path(wd, "Infected.scd.rds"))

rm(scd)

#DPre
Idents(scd.DPre) <- 'Infected'
scd.markers <- FindAllMarkers(scd.DPre,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scd.markers <- subset(scd.markers, p_val_adj < 0.05)
scd.markers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> scd.markers

scd.Heatmap <- DoHeatmap(scd.DPre, features = scd.markers$gene, group.by = 'ident', group.colors = c("red","blue","grey"), raster = F)

file.name <- paste('RNA.DPre.infected','pdf', sep='.')
pdf( file = file.path(wd, 'write.2', file.name ), width = 15, height = 12)
print( scd.Heatmap )
dev.off()

write.csv(scd.markers, file.path(wd, "write.2", "RNA.DPre.infected.markers.csv"))

#DPsel
Idents(scd.DPsel) <- 'Infected'
scd.markers <- FindAllMarkers(scd.DPsel,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scd.markers <- subset(scd.markers, p_val_adj < 0.05)
scd.markers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> scd.markers

scd.Heatmap <- DoHeatmap(scd.DPsel, features = scd.markers$gene, group.by = 'ident', group.colors = c("red","blue","grey"), raster = F)

file.name <- paste('RNA.DPsel.infected','pdf', sep='.')
pdf( file = file.path(wd, 'write.2', file.name ), width = 15, height = 12)
print( scd.Heatmap )
dev.off()

write.csv(scd.markers, file.path(wd, "write.2", "RNA.DPsel.infected.markers.csv"))

#CD4SP
Idents(scd.CD4SP) <- 'Infected'
scd.markers <- FindAllMarkers(scd.CD4SP,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scd.markers <- subset(scd.markers, p_val_adj < 0.05)
scd.markers %>%
  group_by(cluster) %>%   
  arrange(desc(avg_log2FC), .by_group = TRUE) -> scd.markers

scd.Heatmap <- DoHeatmap(scd.CD4SP, features = scd.markers$gene, group.by = 'ident', group.colors = c("red","blue","grey"), raster = F)

file.name <- paste('RNA.CD4SP.infected','pdf', sep='.')
pdf( file = file.path(wd, 'write.2', file.name ), width = 15, height = 12)
print( scd.Heatmap )
dev.off()

write.csv(scd.markers, file.path(wd, "write.2", "RNA.CD4SP.infected.markers.csv"))

gc()