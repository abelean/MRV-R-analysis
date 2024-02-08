#Apply new DN/DP labels to main object
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
options(future.globals.maxSize = 32000 * 1024^2) #32GB

# source
wd <- getwd()
source( file = file.path(wd, 'src', 'shared_r_functions.R' ) )

# Load the Seurat Object 'tsub' ('Single Cell Data') or 'tsub' ('T cell subset')
rds.path <- file.path(wd, 'scd.lin_3.integrated.new.rds')
scd <- readRDS(rds.path)
scd # Sanity check
head( Idents(object = scd) ) # Sanity check

# Load the Seurat Object 'tsub' ('Single Cell Data') or 'tsub' ('T cell subset')
rds.path <- file.path(wd, 'scd.DP.rds')
scd.DP <- readRDS(rds.path)
scd.DP # Sanity check
head( Idents(object = scd.DP) ) # Sanity check

# Load the Seurat Object 'tsub' ('Single Cell Data') or 'tsub' ('T cell subset')
rds.path <- file.path(wd, 'scd.DN.rds')
scd.DN <- readRDS(rds.path)
scd.DN # Sanity check
head( Idents(object = scd.DN) ) # Sanity check

#Get barcodes and metadata from dp and dn files
DP.list <- data.frame('barcodes' = rownames(scd.DP@meta.data), 'lin_X2' = scd.DP$DP)
DN.list <- data.frame('barcodes' = rownames(scd.DN@meta.data), 'lin_X2' = scd.DN$DN)

print(table(DP.list$lin_X2))
print(table(DN.list$lin_X2)) #Sanity check

D.list <- rbind(DP.list, DN.list)

print(table(D.list$lin_X2)) #Sanity Check

#get metadata from whole object and remove DN/DP values
join.list <- data.frame('barcodes' = rownames(scd@meta.data), 'lin_X2' = scd$lin_XX)

print(table(join.list$lin_X2)) #Sanity check

join.list$lin_X2[join.list$lin_X2 %in% c('DN1','DN2','DN3','DN4','Dpbla','Dpre','Dpel')] <- NA

print(table(join.list$lin_X2)) #Sanity check

#Left join D to barcodes
barcodes <- data.frame("barcodes" = rownames(scd@meta.data))
D.barcodes <- left_join(barcodes, D.list, "barcodes")

print(table(D.barcodes$lin_X2)) #Sanity Check

#coalesce to combine lists
join.list$lin_X2 <- coalesce(join.list$lin_X2, D.barcodes$lin_X2)

print(table(join.list$lin_X2)) #Sanity check

#addmetadata
scd <- AddMetaData(scd, join.list$lin_X2, "lin_X2")

#Save obj
saveRDS(scd, file.path(wd, "scd.lin_3.integrated.new.rds"))

gc()