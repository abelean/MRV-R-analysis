# DoMultiBarHeatmap source: https://github.com/satijalab/seurat/issues/2201
#Identify infected and replicating populations, make ORF heatmaps
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
#library( VGAM )
library( patchwork )
library( ggpubr )
library( harmony )

#Increase memory usage
options(future.globals.maxSize = 32000 * 1024^2) #32GB

# Load the Seurat Object 'tsub' ("Single Cell Data") or 'tsub' ("T cell subset")
wd <- getwd()
source( file = file.path(wd, "src", "shared_r_functions.R" ) )

rds.path <- file.path(wd, 'scd.lin_3.integrated.new.rds')
tsub <- readRDS(rds.path)
tsub # Sanity check
head( Idents(object = tsub) ) # Sanity check

#subset out MRV
scd <- subset(tsub, Infection == 'MRV')
#scd <- subset(tsub, lin_X2 == 'Dpbla')
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
#markers and heatmap
DefaultAssay(scd) <- "virus"

my_levels <- c('DN1/2', 'DN3', 'DN4', 'DPbla', 'DPre', 'DPsel', 'ISP', 'CD4_SP', 'CD8_SP', 'NK/ILC', 'Endo/Mes/cTEC/mTEC', 'DC/ETP/Mac', 'RBC')
scd$lin_X2 <- factor(scd$lin_X2, levels = c(my_levels))

Idents(scd) <- 'lin_X2'
scd <- ScaleData(scd)

scd.markers <- FindAllMarkers(scd)

scd.markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) -> scd.markers

scd.markers <- subset(scd.markers, gene %in% orf.list$gene)
scd.markers <- scd.markers[match(orf.list$gene, scd.markers$gene),]
scd.markers <- na.omit(scd.markers)

###subset U values and redo heatmap
orf.sub <- subset(orf.list, kinetics != "U")
table(orf.list$kinetics)
table(orf.sub$kinetics) #Sanity check

#reorder list by kinetics
x <- c("IE","E","E/L","L")
orf.sub %>%
  mutate(kinetics =  factor(kinetics, levels = x)) %>%
  arrange(kinetics) -> orf.sub


scd.sub.markers <- subset(scd.markers, gene %in% orf.sub$gene)
scd.sub.markers <- scd.sub.markers[match(orf.sub$gene, scd.sub.markers$gene),]

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

###Heatmap of inf vs rep
####define infection####
inf <- subset(scd, ORF55 >= 1 & ORF83 >= 1 & ORF103 >= 1)
rep <- subset(scd,
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
inf.list <- data.frame(barcodes = rownames(scd@meta.data))
rep.list <- data.frame(barcodes = rownames(scd@meta.data))
join.list <- data.frame(barcodes = rownames(scd@meta.data))

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

##addmetadata
scd <- AddMetaData(scd, join.list$Infected, "Infected")
table(scd$Infected)#Sanity check

#### ####
###heatmaps comparing infected vs replicating
DefaultAssay(scd) <- "virus"
Idents(scd) <- 'Infected'

scd.markers <- FindAllMarkers(scd)

scd.markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) -> scd.markers

orf.sub <- subset(orf.list, kinetics != "U")
table(orf.list$kinetics)
table(orf.sub$kinetics) #Sanity check

#reorder list by kinetics
x <- c("IE","E","E/L","L")
orf.sub %>%
  mutate(kinetics =  factor(kinetics, levels = x)) %>%
  arrange(kinetics) -> orf.sub

scd.sub.markers <- subset(scd.markers, gene %in% orf.sub$gene)
scd.sub.markers <- scd.sub.markers[match(orf.sub$gene, scd.sub.markers$gene),]

scd.Heatmap <- DoHeatmap(scd, label = F, features = scd.sub.markers$gene, group.by = 'ident', angle = 90, size = 2.5, raster = F) + 
  theme(text = element_text(size = 6))

file.name <- paste('ORFs.infected','pdf', sep='.')
pdf( file = file.path(wd, 'write.2', file.name ), width = 15 )
print( scd.Heatmap )
dev.off()

#Unknown heat map
DefaultAssay(scd) <- "virus"
Idents(scd) <- 'Infected'
#### ####
unknown <- c("ORF1",
             "ORF2",
             "ORF3",
             "ORF4",
             "ORF5",
             "ORF6",
             "ORF9",
             "ORF10",
             "ORF11",
             "ORF12",
             "ORF13",
             "ORF23",
             "ORF35",
             "ORF36",
             "ORF38",
             "ORF43",
             "ORF44",
             "ORF46",
             "ORF52",
             "ORF60",
             "ORF68",
             "ORF70",
             "ORF75",
             "ORF76",
             "ORF81",
             "ORF82",
             "ORF84",
             "ORF87",
             "ORF90",
             "ORF93",
             "ORF94",
             "ORF96",
             "ORF97",
             "ORF100",
             "ORF101",
             "ORF104",
             "ORF105",
             "ORF106",
             "ORF107",
             "ORF108",
             "ORF109",
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
             "ORF128",
             "ORF125a",
             "ORF126a",
             "ORF127a",
             "ORF128a",
             "U25",
             "U32",
             "ORF125b",
             "ORF126b",
             "ORF127b",
             "ORF128b")
#### ####

scd.Heatmap <- DoHeatmap(scd, label = F, features = unknown, group.by = 'ident', angle = 90, size = 2.5, raster = F) + 
  theme(text = element_text(size = 6))

file.name <- paste('Unknown.infected','pdf', sep='.')
pdf( file = file.path(wd, 'write.2', file.name ), width = 15 )
print( scd.Heatmap )
dev.off()

####Add to main object and save
infected.list <- data.frame(barcodes = rownames(scd@meta.data), Infected = scd$Infected)
infected.join <- data.frame(barcodes = rownames(tsub@meta.data))

infected.join <- left_join(infected.join, infected.list, by = "barcodes")

print(table(infected.join$Infected))
print(table(scd$Infected))#Sanity check

###Addmetadata and label uninfected
tsub <- AddMetaData(tsub, infected.join$Infected, "Infected")
print(table(tsub$Infected))#Sanity check

###plot to look at expression levels
plot <- data.frame(table(colSums(tsub@assays$virus@counts)))
plot2 <- plot(plot$Var1, plot$Freq, xlim = range(0:50), ylim = range(0:100), xlab = "colSums of counts", ylab = "Frequency")

file.name <- paste('uninfected.histogram','pdf', sep='.')
pdf( file = file.path(wd, 'write.2', file.name ), width = 15 )
print( plot2 )
dev.off()

n <- length(tsub$Infected)
for(i in 1:n){
  if(colSums(tsub@assays$virus@counts)[i] == 0){ #select value from plot
    tsub$Infected[i] <- "None"
  } else if(colSums(tsub@assays$virus@counts)[i] < 20){ # remove 0 values to make comparison
    tsub$Infected[i] <- "Un-Infected"
    }
}

print(table(tsub$Infected))#Sanity check

#Save obj
saveRDS(tsub, file.path(wd, "Infected.rds"))

#Make percents
dft <- as.data.frame.matrix(table(tsub$lin_X2,tsub$Infected))
dft$All <- dft$Replicating + dft$Infected
write.csv( dft, file.path( wd, 'Infected_2_Table.csv'))

gc()