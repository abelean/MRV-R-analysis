#Create new annotation for whole object based on TEC analysis
#load libaries
library( ggplot2 )
library( cowplot )
library( dplyr ) 
library( tidyr )
library( purrr )
library( tibble )
library( stringr )
library( Seurat )
library( pheatmap )
library( RColorBrewer )
library( VGAM )
library( readr )
library( nichenetr )
library( future )

#Increase memory usage and set plan
options(future.globals.maxSize = 32 * 1024^2) #32GB
#plan("multicore")

#ID project directory
wd <- getwd()
source( file = file.path(wd, "src", "shared_r_functions.R" ) )

# Load the Seurat Object 'scd' ("Single Cell Data")
rds.path <- file.path(wd, 'scd.lin_3.integrated.new.rds')
scd <- readRDS(rds.path)
scd # Sanity check
head( Idents(object = scd) ) # Sanity check

rds.path <- file.path(wd, 'TEC.rds')
scd.tec <- readRDS(rds.path)
scd.tec # Sanity check
head( Idents(object = scd.tec) ) # Sanity check

##Replace labels
#Get barcodes and metadata from TEC files
TEC.list <- data.frame('barcodes' = rownames(scd.tec@meta.data), 'lin_X3' = scd.tec$TEC)

print(table(TEC.list$lin_X3)) #Sanity check

#get metadata from whole object and remove TEC values
join.list <- data.frame('barcodes' = rownames(scd@meta.data), 'lin_X3' = scd$lin_X2)

print(table(join.list$lin_X3)) #Sanity check

join.list$lin_X3[join.list$lin_X3 %in% c('Endo/Mes/cTEC/mTEC', 'DC/ETP/Mac')] <- NA

print(table(join.list$lin_X3)) #Sanity check

#Left join D to barcodes
barcodes <- data.frame("barcodes" = rownames(scd@meta.data))
TEC.barcodes <- left_join(barcodes, TEC.list, "barcodes")

print(table(TEC.barcodes$lin_X3)) #Sanity Check

#coalesce to combine lists
join.list$lin_X3 <- coalesce(join.list$lin_X3, TEC.barcodes$lin_X3)

print(table(join.list$lin_X3)) #Sanity check

#Replace DP/DN
join.list$lin_X3[join.list$lin_X3 %in% c("DN1/2", "DN3", "DN4")] <- "DN"
join.list$lin_X3[join.list$lin_X3 %in% c("DPbla", "DPre", "DPsel")] <- "DP"

print(table(join.list$lin_X3)) #Sanity check
#addmetadata
scd <- AddMetaData(scd, join.list$lin_X3, "lin_X3")

#Save obj
saveRDS(scd, file.path(wd, "scd.lin_3.integrated.new.rds"))

gc()