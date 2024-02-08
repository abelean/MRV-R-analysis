#Graphs of metadata, and some ORF expression
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
library(reshape2)

#Increase memory usage
options(future.globals.maxSize = 32000 * 1024^2) #32GB

# Load the Seurat Object 'tsub' ("Single Cell Data") or 'tsub' ("T cell subset")
wd <- getwd()
rds.path <- file.path(wd, 'scd.lin_3.integrated.new.rds')
tsub <- readRDS(rds.path)
tsub # Sanity check
head( Idents(object = tsub) ) # Sanity check

rds.path <- file.path(wd, 'scd.lin_3.rds')
scd.raw <- readRDS(rds.path)
scd.raw # Sanity check
head( Idents(object = scd.raw) ) # Sanity check

#levels
my_levels <- c('DN1/2', 'DN3', 'DN4', 'ISP', 'DPbla', 'DPre', 'DPsel', 'CD4_SP', 'CD8_SP', 'NK/ILC', 'Endo/Mes/cTEC/mTEC', 'DC/ETP/Mac', 'RBC')
tsub$lin_X2 <- factor(tsub$lin_X2, levels = c(my_levels))

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

###Proportion bar graph
#####Calculations #####change to relative and get cell counts ####print heatmap data
table <- table(tsub$Infection, tsub$lin_X2)
dft <- as.data.frame(rbind(table))
write.csv( dft, file.path( wd, 'MRV.Mock_Table.csv'))

#plot calculations
dft2 <- dft
dft2$Infection <- rownames(dft)
dft2 <- melt(dft2)
plot <- ggplot(dft2, aes(fill = variable, y = value, x = Infection)) +  scale_fill_manual(values = c(colors)) +  geom_bar(stat = "identity", position = "fill") + 
 scale_y_continuous() + theme(panel.border = element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank())

file.name <- paste( 'MRV.Mock_Bar','pdf', sep='.')
pdf( file = file.path(wd, 'write', file.name ) )
print( plot )
dev.off()

###MRV mock UMAP
Idents(tsub) <- "lin_X2"
plot <- UMAPPlot(tsub, split.by = "Infection", cols = colors)

file.name <- paste( 'MRV.Mock_UMAP','pdf', sep='.')
pdf( file = file.path(wd, 'write', file.name ) )
print( plot )
dev.off()

###MRV orf feature plots
tsub.mrv <- subset(tsub, Infection == "MRV")
DefaultAssay(tsub.mrv) <- "virus"
tsub.mrv <- NormalizeData(tsub.mrv)

Idents(tsub.mrv) <- "lin_X2"
my_levels <- c('DN1/2', 'DN3', 'DN4', 'ISP', 'DPbla', 'DPre', 'DPsel', 'CD4_SP', 'CD8_SP', 'NK/ILC', 'Endo/Mes/cTEC/mTEC', 'DC/ETP/Mac', 'RBC')

tsub.mrv@active.ident <- factor(x = tsub.mrv@active.ident, levels = my_levels)

plot <- VlnPlot(tsub.mrv, c("ORF103","ORF83","ORF69"), pt.size = 0, ncol = 2, cols = colors) 

file.name <- paste( 'MRV_ORF_Violin','pdf', sep='.')
pdf( file = file.path(wd, 'write', file.name ))
print( plot )
dev.off()

gc()