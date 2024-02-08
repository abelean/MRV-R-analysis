#Comparing genes identified as unregulated in replicating cells
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
library( gplots )
library( enrichR )
library( reshape2 ) 

#Increase memory usage
options(future.globals.maxSize = 32000 * 1024^2) #32GB

# Load the Seurat Object 'scd' ("Single Cell Data") or 'scd' ("T cell subset")
wd <- getwd()
source( file = file.path(wd, "src", "shared_r_functions.R" ) )

rds.path <- file.path(wd, 'scd.lin_3.integrated.new.rds')
scd <- readRDS(rds.path)
scd # Sanity check
head( Idents(object = scd) ) # Sanity check


CD4SP.markers <- read.csv(file.path(wd, "write.2", "RNA.CD4SP.infected.markers.csv"))
DPre.markers <- read.csv(file.path(wd, "write.2", "RNA.DPre.infected.markers.csv"))
DPsel.markers <- read.csv(file.path(wd, "write.2", "RNA.DPsel.infected.markers.csv"))

CD4SP.markers$X <- NULL
DPre.markers$X <- NULL
DPsel.markers$X <- NULL

#Replicating
CD4SP.rep <- subset(CD4SP.markers, cluster == "Replicating")
DPre.rep <- subset(DPre.markers, cluster == "Replicating")
DPsel.rep <- subset(DPsel.markers, cluster == "Replicating")

rep <- list(CD4SP = CD4SP.rep$gene,
            DPre = DPre.rep$gene,
            DPsel = DPsel.rep$gene)

file.name <- paste('Venn.Replicating', 'pdf', sep='.')
pdf( file = file.path(wd, 'write.2', file.name ))
print( venn(rep) )
dev.off()

table <- intersect(intersect(CD4SP.rep$gene, DPre.rep$gene), DPsel.rep$gene)
write.csv(table, file.path(wd, "write.2", "Venn.Intersect.Replicating.csv"))
print(length(intersect(intersect(CD4SP.rep$gene, DPre.rep$gene), DPsel.rep$gene)))#Sanity check


###General Enrichr pipeline
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes

websiteLive <- TRUE

dbs <- listEnrichrDbs()

if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

dbs <- c("WikiPathways_2019_Mouse", "HDSigDB_Mouse_2021", "KEGG_2019_Mouse")

csv.list <- list()

#Function for wait
wait <- function(x)
{
  p1 <- proc.time()
  Sys.sleep(x)
  proc.time() - p1 # The cpu usage should be negligible
}

###make pathway plots
wait(10)
enriched <- c()
if (websiteLive) {
  enriched <- enrichr(c(table), dbs) #copy genes from heatmaps here
}

if (websiteLive) enriched[["WikiPathways_2019_Mouse"]]

enriched[[1]] <- enriched[[1]] %>% 
  mutate(num = str_count(Genes, ";")+1)

enriched[[1]] <- filter(enriched[[1]], Adjusted.P.value < 0.05)
enriched[[1]] %>% 
  arrange(Adjusted.P.value) -> enriched[[1]]

file.name <- paste( "Replicating.Intersect", 'csv', sep='.')
write.csv(enriched[[1]], file.path(wd, 'write.2', 'Raw', file.name))

names(enriched)[1] <- "Replicating.Intersect"
csv.list[[1]] <- enriched[1]
enrich <- as.data.frame(csv.list[1])
colnames(enrich) <- colnames(csv.list[[1]][[1]])
enrich$Term  <- gsub(" WP.*","", enrich$Term)
enrich %>% 
  arrange(desc(Adjusted.P.value)) -> enrich
enrich$Term <- factor(enrich$Term, levels = enrich$Term)

plot <- ggplot(enrich, aes(x = Term, y = num)) +
  geom_col(aes(fill = Adjusted.P.value), width = 0.7) + 
  coord_flip() +
  ggtitle(names(enriched[1])) +
  ylab("Genes") +
  xlab("Pathway") +
  scale_fill_gradient(low="red",high="blue", trans = 'log', breaks=c(0.00005,0.0005,0.005,0.05),limits=c(min(enrich$Adjusted.P.value),0.051)) +
  theme_classic()

file.name <- paste( "Replicating.Intersect", 'pdf', sep='.')
pdf( file = file.path(wd, 'write.2', 'EnrichR', file.name ) )
print( plot )
dev.off() 

## 2 intersects+
list1 <- intersect(CD4SP.rep$gene, DPre.rep$gene)
list2 <- intersect(CD4SP.rep$gene, DPsel.rep$gene)
list3 <- intersect(DPre.rep$gene, DPsel.rep$gene)

table2 <- c(list1, list2, list3)
table2 <- unique(table2)

print(117+35+76+43)
print(length(table2))

write.csv(table2, file.path(wd, "write.2", "Venn.Intersect.All.Replicating.csv"))

###make pathway plots
wait(10)
enriched <- c()
if (websiteLive) {
  enriched <- enrichr(c(table2), dbs) #copy genes from heatmaps here
}

if (websiteLive) enriched[["WikiPathways_2019_Mouse"]]

enriched[[1]] <- enriched[[1]] %>% 
  mutate(num = str_count(Genes, ";")+1)

enriched[[1]] <- filter(enriched[[1]], Adjusted.P.value < 0.05)
enriched[[1]] %>% 
  arrange(Adjusted.P.value) -> enriched[[1]]

file.name <- paste( "Replicating.Intersect.All", 'csv', sep='.')
write.csv(enriched[[1]], file.path(wd, 'write.2', 'Raw', file.name))

names(enriched)[1] <- "Replicating.Intersect.All"
csv.list[[1]] <- enriched[1]
enrich <- as.data.frame(csv.list[1])
colnames(enrich) <- colnames(csv.list[[1]][[1]])
enrich$Term  <- gsub(" WP.*","", enrich$Term)
enrich %>% 
  arrange(desc(Adjusted.P.value)) -> enrich
enrich$Term <- factor(enrich$Term, levels = enrich$Term)

plot <- ggplot(enrich, aes(x = Term, y = num)) +
  geom_col(aes(fill = Adjusted.P.value), width = 0.7) + 
  coord_flip() +
  ggtitle(names(enriched[1])) +
  ylab("Genes") +
  xlab("Pathway") +
  scale_fill_gradient(low="red",high="blue", trans = 'log', breaks=c(0.00005,0.0005,0.005,0.05),limits=c(min(enrich$Adjusted.P.value),0.051)) +
  theme_classic() +
  scale_y_continuous(limits = c(0,9), breaks = c(0,3,6,9))

file.name <- paste( "Replicating.Intersect.All", 'pdf', sep='.')
pdf( file = file.path(wd, 'write.2', 'EnrichR', file.name ) )
print( plot )
dev.off()

##Dotplot
# Load the Seurat Object 'scd' ("Single Cell Data")
rds.path <- file.path(wd, 'Infected.scd.rds')
scd.inf <- readRDS(rds.path)
scd.inf # Sanity check
head( Idents(object = scd.inf) ) # Sanity check

cholesterol <- unique(c("Rpe","Pgd",
                        "Fdps","Idi1","Hmgcs1","Msmo1","Hmgcr","Cyp51",
  "Idi1","Fdps","Hmgcs1","Msmo1","Hmgcr","Fads1","Srebf2","Cyp51"))
cell.cycle <- unique(c("Prim2","Cdc45","Rfc4","Rfc2","Rpa3","Gmnn","Rpa2","Mcm5","Mcm6",
  "Prim2","Cdc45","Rpa3","E2f1","Cdk1","Rpa2","Mcm5","Mcm6",
  "Dhfr","Rrm1","Rrm2"))


genes <- c(cell.cycle, cholesterol)

scd.dot <- subset(scd.inf, Infection == "MRV")
scd.dot <- subset(scd.dot, lin_X2 %in% c("DPre","DPsel","CD4_SP"))

DefaultAssay(scd.dot) <- "RNA"
scd.dot$Infected <- factor(scd.dot$Infected, c("Un-Infected", "Infected", "Replicating"))
Idents(scd.dot) <- "Infected"
plot <- DotPlot(scd.dot, features = rev(genes), cols = "RdYlBu", split.by = "lin_X2") + coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

file.name <- paste('Dotplot.Replicating','pdf', sep='.')
pdf( file = file.path( wd , 'write.2', file.name ) , width = 10 , height = 10 )
print(plot)
dev.off()

gc()