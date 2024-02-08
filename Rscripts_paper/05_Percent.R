#Calculatin frequencies of infected/replicating in MRV cell types
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
#library( monocle3 )
library( VGAM )
library( patchwork )
library( ggpubr )

#Increase memory usage
options(future.globals.maxSize = 32000 * 1024^2) #32GB

# Load the Seurat Object 'tsub' ("Single Cell Data") or 'tsub' ("T cell subset")
wd <- getwd()
rds.path <- file.path(wd, 'scd.lin_3.integrated.new.rds')
scd.raw <- readRDS(rds.path)
scd.raw # Sanity check
head( Idents(object = scd.raw) ) # Sanity check

#Levels
my_levels <- c('DN1/2', 'DN3', 'DN4', 'ISP', 'DPbla', 'DPre', 'DPsel', 'CD4_SP', 'CD8_SP', 'NK/ILC', 'Endo/Mes/cTEC/mTEC', 'DC/ETP/Mac', 'RBC')
scd.raw$lin_X2 <- factor(scd.raw$lin_X2, levels = c(my_levels))

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

#subset by presence of genes to split into infections and replication
DefaultAssay(scd.raw) <- 'virus'
inf <- subset(scd.raw, ORF55 >= 1 & ORF83 >= 1 & ORF103 >= 1)
rep <- subset(scd.raw,
              ( ORF53>=3 & ORF69>=3 & ORF73>=3 & ORF86>=3 ) |
                ( ORF45>=3 & ORF69>=3 & ORF73>=3 & ORF86>=3 ) |
                ( ORF45>=3 & ORF53>=3 & ORF73>=3 & ORF86>=3 ) |
                ( ORF45>=3 & ORF53>=3 & ORF69>=3 & ORF86>=3 ) |
                ( ORF45>=3 & ORF53>=3 & ORF69>=3 & ORF73>=3 ))

###calculate percentages for MRV vs Mock
##load file
rds.path <- file.path(wd, 'Infected.rds')
scd.raw <- readRDS(rds.path)
scd.raw # Sanity check
head( Idents(object = scd.raw) ) # Sanity check

##Subset objects
Idents(scd.raw) <- "Infected"
inf <- subset(scd.raw, Infected %in% 'Infected')
rep <- subset(scd.raw, Infected %in% c('Both','Replicating'))

#Calculate percentages for MRV
scd.mrv <- subset(scd.raw, Infection == "MRV")

#Calculate percentages
#Create dataframe to input calculations into
x <- as.data.frame(table(scd.mrv$lin_X2))
head(x)
table(scd.mrv$lin_X2)#sanity check

x <- cbind(x, 'infection')
x <- cbind(x, 'replication')

#Create data for infection and replication, and add empty rows where necessary using a left join
ii <- as.data.frame(table(inf$lin_X2))
dfi <- as.data.frame(seq(1, 13))
colnames(dfi) <- 'Var1'
dfi$Var1 <- x$Var1
dfi <- left_join(dfi, ii, by = 'Var1')
dfi[is.na(dfi)] <- 0
head(c(ii, dfi)) # Sanity check

r <- as.data.frame(table(rep$lin_X2))
dfr <- as.data.frame(seq(1, 13))
colnames(dfr) <- 'Var1'
dfr$Var1 <- x$Var1
dfr <- left_join(dfr, r, by = 'Var1')
dfr[is.na(dfr)] <- 0
head(c(r, dfr)) # Sanity check

#Do calculations
for(i in 1:length(x$Freq)){
  x$`"infection"`[i] <- (dfi$Freq[i] / x$Freq[i]) * 100
  x$`"replication"`[i] <- (dfr$Freq[i] / x$Freq[i]) * 100
}
#add freq for inf and rep
x$`"freq_inf"` <- dfi$Freq
x$`"freq_rep"` <- dfr$Freq
x <- x[match(my_levels, x$Var1),]
rownames(x) <- x$Var1
x$Var1 <- NULL

#Save calculations(calculate in excel)
write.csv( x , file.path( wd, 'Percent.MRV.csv' ))

gc()