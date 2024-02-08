# Monocle3 workflow based on https://ucdavis-bioinformatics-training.github.io/2021-August-Advanced-Topics-in-Single-Cell-RNA-Seq-Trajectory-and-Velocity/data_analysis/monocle_fixed
#Running pseudo time on infected/replicating populations, then generating pseudo kinetics plots
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
library( SeuratWrappers )
library( monocle3 )
#library( VGAM )
library( patchwork )
library( ggpubr )

#Increase memory usage
options(future.globals.maxSize = 32000 * 1024^2) #32GB

# Load the Seurat Object 'scd' ("Single Cell Data")
wd <- getwd()
rds.path <- file.path(wd, 'Infected.rds')
scd.raw1 <- readRDS(rds.path)
scd.raw1 # Sanity check
head( Idents(object = scd.raw1) ) # Sanity check

#### Redo dim reduction based on viral####
DefaultAssay(scd.raw1) <- "virus"

scd.raw <- subset(scd.raw1, Infected != "None")
scd.raw <- subset(scd.raw, Infected != "Un-Infected")
scd.raw <- FindVariableFeatures(scd.raw, assay = "virus")

scd.raw <- NormalizeData(scd.raw, assay = "virus")
all.genes <- rownames(scd.raw)
scd.raw<- ScaleData(scd.raw, features = all.genes, assay = "virus")

#<- FindNeighbors(scd.raw, assay = "virus", reduction = "pca")

scd.raw<- ScaleData(scd.raw, features = all.genes, assay = "virus")
#Run PCA on filtered data
scd.raw <- RunPCA(scd.raw, features = VariableFeatures(scd.raw), assay = "virus")
#Run Harmony on processed data

#scd <- RunHarmony(scd.raw, group.by.vars = "lin_XX", assay.use = "virus") ###group.by???
#Run Umap
scd <- RunUMAP(scd.raw, reduction = "pca", dims = 1:30, verbose = FALSE, assay = "virus") 

Idents(scd) <- "Infected"
plot <- UMAPPlot(scd, cols = c("red","blue")) + theme(legend.text=element_text(size=20))


file.name <- 'Infected.UMAP.pdf'
pdf(file = file.path(wd, "write.5",file.name))
print(plot)
dev.off()

###Add new metadata
meta <- scd$lin_X2
meta[meta %in% c("DN1/2","DN3","DN4")] <- "DN"
meta[meta %in% c("DPbla","DPre","DPsel")] <- "DP"
meta[meta %in% c("DC/ETP/Mac","Endo/Mes/cTEC/mTEC","ISP","RBC")] <- "Other"

scd <- AddMetaData(scd, meta, "lin_X3")

print(table(scd$lin_X2))
print(table(scd$lin_X3))

# Convert Seurat Objects to CellDataSet (cds) for Monocle3

cds <- as.cell_data_set(scd)

###estimate_size factor to fix expression ,add gene rownames
cds <- cds[,Matrix::colSums(exprs(cds)) != 0]
cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(scd[["virus"]])

cds <- cds[,Matrix::colSums(exprs(cds)) != 0]
cds <- cluster_cells(cds, resolution=1e-3)#Find resolution that work best

p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)
#Subsetting partions
#integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
#cds <- as.cell_data_set(integrated.sub)
#Trajectory analysis
cds <- learn_graph(cds, use_partition = FALSE, verbose = FALSE)
plot <- plot_cells(cds,
                   color_cells_by = "cluster",
                   label_groups_by_cluster=FALSE,
                   label_leaves=FALSE,
                   label_branch_points=FALSE,
                   label_principal_points = TRUE)#Find root here

file.name <- 'Infected.Root.pdf'
pdf(file = file.path(wd, "write.5",file.name))
print(plot)
dev.off()

# Pseudotime
cds <- order_cells(cds, root_pr_nodes = "Y_62")#select root here, change setting later
plot <- plot_cells(cds,
                   color_cells_by = "pseudotime",
                   group_cells_by = "partition",
                   label_cell_groups = FALSE,
                   label_groups_by_cluster=FALSE,
                   label_leaves=FALSE,
                   label_branch_points=FALSE,
                   label_roots = FALSE,
                   trajectory_graph_color = "grey60")

file.name <- 'Infected.Pseudotime.pdf'
pdf(file = file.path(wd, "write.5",file.name))
print(plot)
dev.off()

##plot genes by pseudotime
#EI
ORF_genes <- c("ORF30", "ORF31", "ORF103")
ORF_lineage_cds <- cds[rowData(cds)$gene_short_name %in% ORF_genes]

plot <-  plot_genes_in_pseudotime(ORF_lineage_cds, color_cells_by = "lin_X3", min_expr = .1, ncol = 3)  + scale_color_brewer(palette = "Set1")

file.name <- 'Pseudotime.Early.Intermediate.pdf'
pdf(file = file.path(wd, "write.5",file.name), height = 2.5)
print(plot)
dev.off()

#E
ORF_genes <- c("ORF21", "ORF26", "ORF29", "ORF41", "ORF47", "ORF50", "ORF51", "ORF53", "ORF54", "ORF55", "ORF57", "ORF64", "ORF65", "ORF66", "ORF72", "ORF74", "ORF88", "ORF91", "ORF98", "ORF102", "ORF110")
ORF_lineage_cds <- cds[rowData(cds)$gene_short_name %in% ORF_genes]

plot <-  plot_genes_in_pseudotime(ORF_lineage_cds, color_cells_by = "lin_X3", min_expr = .1, ncol = 3) + scale_color_brewer(palette = "Set1")

file.name <- 'Pseudotime.Early.pdf'
pdf(file = file.path(wd, "write.5",file.name))
print(plot)
dev.off()

#E/L
ORF_genes <- c("ORF7",
               "ORF8",
               "ORF14",
               "ORF15",
               "ORF16",
               "ORF17",
               "ORF18",
               "ORF19",
               "ORF27",
               "ORF39",
               "ORF40",
               "ORF45",
               "ORF56",
               "ORF58",
               "ORF59",
               "ORF62",
               "ORF71",
               "ORF73",
               "ORF80",
               "ORF83",
               "ORF85",
               "ORF93")
ORF_lineage_cds <- cds[rowData(cds)$gene_short_name %in% ORF_genes]

plot <-  plot_genes_in_pseudotime(ORF_lineage_cds, color_cells_by = "lin_X3", min_expr = .1, ncol = 3) + scale_color_brewer(palette = "Set1")

file.name <- 'Pseudotime.Early.Late.pdf'
pdf(file = file.path(wd, "write.5",file.name))
print(plot)
dev.off()

#L
ORF_genes <- c("ORF20",
               "ORF22",
               "ORF24",
               "ORF25",
               "ORF28",
               "ORF33",
               "ORF34",
               "ORF37",
               "ORF42",
               "ORF48",
               "ORF49",
               "ORF61",
               "ORF63",
               "ORF67",
               "ORF69",
               "ORF77",
               "ORF78",
               "ORF79",
               "ORF86",
               "ORF89",
               "ORF92",
               "ORF99")
ORF_lineage_cds <- cds[rowData(cds)$gene_short_name %in% ORF_genes]

plot <-  plot_genes_in_pseudotime(ORF_lineage_cds, color_cells_by = "lin_X3", min_expr = .1, ncol = 3) + scale_color_brewer(palette = "Set1")

file.name <- 'Pseudotime.Late.pdf'
pdf(file = file.path(wd, "write.5",file.name))
print(plot)
dev.off()

#U1
ORF_genes <- c("ORF1",
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
               "ORF93")
ORF_lineage_cds <- cds[rowData(cds)$gene_short_name %in% ORF_genes]

plot <-  plot_genes_in_pseudotime(ORF_lineage_cds, color_cells_by = "lin_X3", min_expr = .1, ncol = 3) + scale_color_brewer(palette = "Set1")

file.name <- 'Pseudotime.Unknown1.pdf'
pdf(file = file.path(wd, "write.5",file.name))
print(plot)
dev.off()

#U2
ORF_genes <- c("ORF94",
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
               "ORF128")
ORF_lineage_cds <- cds[rowData(cds)$gene_short_name %in% ORF_genes]

plot <-  plot_genes_in_pseudotime(ORF_lineage_cds, color_cells_by = "lin_X3", min_expr = .1, ncol = 3) + scale_color_brewer(palette = "Set1")

file.name <- 'Pseudotime.Unknown2.pdf'
pdf(file = file.path(wd, "write.5",file.name))
print(plot)
dev.off()

#U3
ORF_genes <- c("ORF125a", "ORF126a", "ORF127a", "ORF128a", "U25", "U32", "ORF125b", "ORF126b", "ORF127b", "ORF128b")
ORF_lineage_cds <- cds[rowData(cds)$gene_short_name %in% ORF_genes]

plot <-  plot_genes_in_pseudotime(ORF_lineage_cds, color_cells_by = "lin_X3", min_expr = .1, ncol = 3) + scale_color_brewer(palette = "Set1")

file.name <- 'Pseudotime.Unknown3.pdf'
pdf(file = file.path(wd, "write.5",file.name))
print(plot)
dev.off()

####Representatives by clustering
#IE
ORF_genes <- c("ORF103","ORF31")
ORF_lineage_cds <- cds[rowData(cds)$gene_short_name %in% ORF_genes]

plot <-  plot_genes_in_pseudotime(ORF_lineage_cds, color_cells_by = "lin_X3", min_expr = .1, ncol = 3) + scale_color_brewer(palette = "Set1")

file.name <- 'Rep.IE.pdf'
pdf(file = file.path(wd, "write.5",file.name), height = 2.5)
print(plot)
dev.off()

#E
ORF_genes <- c("ORF83","ORF56")
ORF_lineage_cds <- cds[rowData(cds)$gene_short_name %in% ORF_genes]

plot <-  plot_genes_in_pseudotime(ORF_lineage_cds, color_cells_by = "lin_X3", min_expr = .1, ncol = 3) + scale_color_brewer(palette = "Set1")

file.name <- 'Rep.E.pdf'
pdf(file = file.path(wd, "write.5",file.name), height = 2.5)
print(plot)
dev.off()

#EL
ORF_genes <- c("ORF73","ORF79")
ORF_lineage_cds <- cds[rowData(cds)$gene_short_name %in% ORF_genes]

plot <-  plot_genes_in_pseudotime(ORF_lineage_cds, color_cells_by = "lin_X3", min_expr = .1, ncol = 3) + scale_color_brewer(palette = "Set1")

file.name <- 'Rep.EL.pdf'
pdf(file = file.path(wd, "write.5",file.name), height = 2.5)
print(plot)
dev.off()

#L
ORF_genes <- c("ORF69","ORF72")
ORF_lineage_cds <- cds[rowData(cds)$gene_short_name %in% ORF_genes]

plot <-  plot_genes_in_pseudotime(ORF_lineage_cds, color_cells_by = "lin_X3", min_expr = .1, ncol = 3) + scale_color_brewer(palette = "Set1")

file.name <- 'Rep.L.pdf'
pdf(file = file.path(wd, "write.5",file.name), height = 2.5)
print(plot)
dev.off()

###Infected
#IE
ORF_genes <- c("ORF103","ORF31")
ORF_lineage_cds <- cds[rowData(cds)$gene_short_name %in% ORF_genes]

plot <-  plot_genes_in_pseudotime(ORF_lineage_cds, color_cells_by = "Infected", min_expr = .1, ncol = 3) + scale_color_brewer(palette = "Set1")

file.name <- 'Rep.IE.Infected.pdf'
pdf(file = file.path(wd, "write.5",file.name), height = 2.5)
print(plot)
dev.off()

#E
ORF_genes <- c("ORF83","ORF56")
ORF_lineage_cds <- cds[rowData(cds)$gene_short_name %in% ORF_genes]

plot <-  plot_genes_in_pseudotime(ORF_lineage_cds, color_cells_by = "Infected", min_expr = .1, ncol = 3) + scale_color_brewer(palette = "Set1")

file.name <- 'Rep.E.Infected.pdf'
pdf(file = file.path(wd, "write.5",file.name), height = 2.5)
print(plot)
dev.off()

#EL
ORF_genes <- c("ORF73","ORF79")
ORF_lineage_cds <- cds[rowData(cds)$gene_short_name %in% ORF_genes]

plot <-  plot_genes_in_pseudotime(ORF_lineage_cds, color_cells_by = "Infected", min_expr = .1, ncol = 3) + scale_color_brewer(palette = "Set1")

file.name <- 'Rep.EL.Infected.pdf'
pdf(file = file.path(wd, "write.5",file.name), height = 2.5)
print(plot)
dev.off()

#L
ORF_genes <- c("ORF69","ORF72")
ORF_lineage_cds <- cds[rowData(cds)$gene_short_name %in% ORF_genes]

plot <-  plot_genes_in_pseudotime(ORF_lineage_cds, color_cells_by = "Infected", min_expr = .1, ncol = 3) + scale_color_brewer(palette = "Set1")

file.name <- 'Rep.L.Infected.pdf'
pdf(file = file.path(wd, "write.5",file.name), height = 2.5)
print(plot)
dev.off()

gc()