#Generating graphs of expression in specific genes
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

# Load the Seurat Object 'tsub' ("Single Cell Data") or 'tsub' ("T cell subset")
wd <- getwd()
rds.path <- file.path(wd, 'scd.lin_3.integrated.new.rds')
scd <- readRDS(rds.path)
scd # Sanity check
head( Idents(object = scd) ) # Sanity check

Idents(scd) <- "lin_X2"
####Mock####
scd.mock <- subset(scd, Infection == "Mock")

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

scd.DP <- subset(scd.mock, lin_X2 == c("DPbla","DPre","DPsel"))

Idents(scd.DP) <- "lin_X2"
plot <- VlnPlot(scd.DP, c("Cd2", "Mki67"), cols = colors) & ylim(-1,5)

file.name <- 'Violin.DP.Ki67.vs.Cd2.pdf'
pdf(file = file.path(wd, "write.8",file.name))
print(plot)
dev.off()

####MRV####
scd.mrv <- subset(scd, Infection == "MRV")

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

scd.DP <- subset(scd.mrv, lin_X2 == c("DPbla","DPre","DPsel"))

Idents(scd.DP) <- "lin_X2"
plot <- VlnPlot(scd.DP, c("Cd2", "Mki67"), cols = colors) & ylim(-1,5)

file.name <- 'Violin.DP.MRV.Ki67.vs.Cd2.pdf'
pdf(file = file.path(wd, "write.8",file.name))
print(plot)
dev.off()

#### dotplots of cyto ####
scd.cyto <- subset(scd, Infection == "MRV")

scd.cyto$lin_X3 <- factor(scd.cyto$lin_X3, c("ETP","DN","ISP","DP","CD4_SP","CD8_SP","NK/ILC",
                                             "cTEC","mTEC",'DC',"Mono/Mac","Granulocyte","Endo/Mes","RBC"))
Idents(scd.cyto) <- "lin_X3"
DefaultAssay(scd.cyto) <- "RNA"

plot <- DotPlot(scd.cyto, features = rev(c("Ifna1","Il36a","Ifnb","Il7","Il18","Ifng","Il4","Il15","Np","Il12a","Il12b","Tnfsf15","Il2")), 
                cols = "RdYlBu",
                scale.min = 0, scale.max = 100) + coord_flip() +
  theme(axis.text = element_text(size = 15,angle = 45, hjust = 1), legend.position = "right")

file.name <- 'Cytokines.pdf'
pdf(file = file.path(wd, "write.8",file.name))
print(plot)
dev.off()

gc() 