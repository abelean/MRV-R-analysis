# Paley Neurosarcoid single-cell RNA-Seq 
# Shared R functions

#########################################################
# Setup project paths and make directories if necessary #
#########################################################

# path setup
library( here )
library( ggplot2 )
library( cowplot )
library( dplyr ) 
library( tidyr )
library( purrr )
library( tibble )
library( stringr )
library( readr )
library( RColorBrewer )
library( grid )

# paths
figure_path <- file.path(wd, "results", "figures" )
object_path <- file.path(wd, "results", "objects" )
doublet_path <- file.path(wd, "results", "doublets" )
doublet_fig_path <- file.path(wd, "results", "figures", "doublets" )
annotation_path <- file.path(wd, 'annotation' )

# colors for graphs
tissue.colors <- c( Aqueous = 'blue' , Blood = 'red')
T.colors <- c( brewer.pal( 8, "Dark2" ) , brewer.pal( 8, "Pastel2" ) , brewer.pal( 6, "Set2" )  )
subject.colors <- c( brewer.pal( 12, "Paired" ) , brewer.pal( 12, "Set3" ) , brewer.pal( 4, "Accent" )  )
CD4T.colors <- c( brewer.pal( 8, "Dark2" ) , brewer.pal( 8, "Pastel2" )[1:4] )
CD8T.colors <- c( brewer.pal( 8, "Pastel2" )[5:8] , brewer.pal( 6, "Set2" ) )

# setup directories
if (dir.exists(file.path(wd, 'results') ) ) {
    cat("results exists and is a directory.\n")
} else {
    cat("results does not exist - creating.\n")
    dir.create( file.path(wd, 'results') )
}
if (dir.exists(file.path(figure_path) ) ) {
    cat("figures exists and is a directory.\n")
} else {
    cat("figures does not exist - creating.\n")
    dir.create( file.path(figure_path) )
}
if (dir.exists(file.path(object_path) ) ) {
    cat("objects exists and is a directory.\n")
} else {
    cat("objects does not exist - creating.\n")
    dir.create( file.path(object_path) )
}
if (dir.exists(file.path(doublet_path) ) ) {
    cat("doublets exists and is a directory.\n")
} else {
    cat("doublets does not exist - creating.\n")
    dir.create( file.path(doublet_path) )
}
if (dir.exists(file.path(doublet_fig_path) ) ) {
    cat("figures/doublets exists and is a directory.\n")
} else {
    cat("figures/doublets does not exist - creating.\n")
    dir.create( file.path(doublet_fig_path) )
}
if (dir.exists(file.path(annotation_path) ) ) {
    cat("annotation exists and is a directory.\n")
} else {
    cat("annotation does not exist - creating.\n")
    dir.create( file.path(annotation_path) )
}

if (dir.exists(file.path(wd, 'results','calculations') ) ) {
  cat("calculations exists and is a directory.\n")
} else {
  cat("calculations does not exist - creating.\n")
  dir.create( file.path(wd, 'results','calculations') )
}
if (dir.exists(file.path(wd, 'results','calculations','VDJ') ) ) {
  cat("calculations/VDJ exists and is a directory.\n")
} else {
  cat("calculations/VDJ does not exist - creating.\n")
  dir.create( file.path(wd, 'results','calculations','VDJ') )
}
if (dir.exists(file.path(wd, 'results', 'figures') ) ) {
  cat("results/figures exists and is a directory.\n")
} else {
  cat("results/figures does not exist - creating.\n")
  dir.create( file.path(wd, 'results', 'figures') )
}
if (dir.exists(file.path(wd, 'results', 'figures', 'vdj_processing') ) ) {
  cat("results/figures/vdj_processing exists and is a directory.\n")
} else {
  cat("results/figures/vdj_processing does not exist - creating.\n")
  dir.create( file.path(wd, 'results', 'figures', 'vdj_processing') )
}
if (dir.exists(file.path(wd, 'results', 'bcr_qc') ) ) {
  cat("results/bcr_qc exists and is a directory.\n")
} else {
  cat("results/bcr_qc does not exist - creating.\n")
  dir.create( file.path(wd, 'results', 'bcr_qc') )
}
if (dir.exists(file.path(wd, 'results', 'tcr_qc') ) ) {
    cat("results/tcr_qc exists and is a directory.\n")
} else {
    cat("results/tcr_qc does not exist - creating.\n")
    dir.create( file.path(wd, 'results', 'tcr_qc') )
}
if (dir.exists(file.path(wd, 'results', 'vdj_meta_data') ) ) {
    cat("results/vdj_meta_data exists and is a directory.\n")
} else {
    cat("results/vdj_meta_data does not exist - creating.\n")
    dir.create( file.path(wd, 'results', 'vdj_meta_data') )
}
if (dir.exists(file.path(wd, 'qc') ) ) {
    cat("qc exists and is a directory.\n")
} else {
    cat("qc does not exist - creating.\n")
    dir.create( file.path(wd, 'qc') )
}
if (dir.exists(file.path(wd, 'results', 'figures', 'UMAP') ) ) {
  cat("UMAP exists and is a directory.\n")
} else {
  cat("UMAP does not exist - creating.\n")
  dir.create( file.path(wd, 'results', 'figures', 'UMAP') )
}
if (dir.exists(file.path(wd, 'results', 'figures', 'Violin') ) ) {
  cat("Violin exists and is a directory.\n")
} else {
  cat("Violin does not exist - creating.\n")
  dir.create( file.path(wd, 'results', 'figures', 'Violin') )
}
if (dir.exists(file.path(wd, 'sct_steps') ) ) {
    cat("sct_steps exists and is a directory.\n")
} else {
    cat("sct_steps does not exist - creating.\n")
    dir.create( file.path(wd, 'sct_steps') )
}
if (dir.exists(file.path(wd, 'results', 'figures', 'variable_clusters') ) ) {
    cat("variable_clusters exists and is a directory.\n")
} else {
    cat("variable_clusters does not exist - creating.\n")
    dir.create( file.path(wd, 'results', 'figures', 'variable_clusters') )
}
if (dir.exists(file.path(wd, 'results', 'figures', 'chosen_clusters') ) ) {
    cat("chosen_clusters exists and is a directory.\n")
} else {
    cat("chosen_clusters does not exist - creating.\n")
    dir.create( file.path(wd, 'results', 'figures', 'chosen_clusters') )
}
if (dir.exists(file.path(wd, 'results', 'DEGs') ) ) {
  cat("DEGs exists and is a directory.\n")
} else {
  cat("DEGs does not exist - creating.\n")
  dir.create( file.path(wd, 'results', 'DEGs') )
}
if (dir.exists(file.path(wd, 'results', 'Monocle2') ) ) {
  cat("Monocle2 exists and is a directory.\n")
} else {
  cat("Monocle2 does not exist - creating.\n")
  dir.create( file.path(wd, 'results', 'Monocle2') )
}
if (dir.exists(file.path(wd, 'results', 'Monocle2', 'CDS') ) ) {
  cat("CDS exists and is a directory.\n")
} else {
  cat("CDS does not exist - creating.\n")
  dir.create( file.path(wd, 'results', 'Monocle2', 'CDS') )
}
if (dir.exists(file.path(wd, 'results', 'Monocle2', 'Disp_Table') ) ) {
  cat("Disp_Table exists and is a directory.\n")
} else {
  cat("Disp_Table does not exist - creating.\n")
  dir.create( file.path(wd, 'results', 'Monocle2', 'Disp_Table') )
}
if ( dir.exists( file.path( wd, 'results', 'Monocle2', 'Disp_Table', 'lin_markers' ) ) ) {
  cat("lin_markers exists and is a directory.\n")
} else {
  cat("lin_markers does not exist - creating.\n")
  dir.create( file.path( wd, 'results', 'Monocle2', 'Disp_Table', 'lin_markers' ) )
}

if ( dir.exists( file.path( wd, 'results', 'figures', 'demultiplex' ) ) ) {
  cat("figures/demultiplex exists and is a directory.\n")
} else {
  cat("figures/demultiplex does not exist - creating.\n")
  dir.create( file.path( wd, 'results', 'figures', 'demultiplex' ) )
}
if ( dir.exists( file.path( wd, 'results', 'calculations', 'demultiplex' ) ) ) {
  cat("calculations/demultiplex exists and is a directory.\n")
} else {
  cat("calculations/demultiplex does not exist - creating.\n")
  dir.create( file.path( wd, 'results', 'calculations', 'demultiplex' ) )
}
if (dir.exists(file.path(wd, 'results', 'TCR_Analysis') ) ) {
  cat("TCR_Analysis exists and is a directory.\n")
} else {
  cat("TCR_Analysis does not exist - creating.\n")
  dir.create( file.path(wd, 'results', 'TCR_Analysis') )
}

if (dir.exists(file.path(wd, 'results', 'figures', 'T.NK') ) ) {
  cat("figures/T.NK exists and is a directory.\n")
} else {
  cat("figures/T.NK does not exist - creating.\n")
  dir.create( file.path(wd, 'results', 'figures', 'T.NK') )
}
if (dir.exists(file.path(wd, 'results', 'figures', 'T.NK', 'variable_clusters') ) ) {
  cat("figures/T.NK/variable_clusters exists and is a directory.\n")
} else {
  cat("figures/T.NK/variable_clusters does not exist - creating.\n")
  dir.create( file.path(wd, 'results', 'figures', 'T.NK', 'variable_clusters') )
}
if (dir.exists(file.path(wd, 'results', 'figures', 'T.NK', 'chosen_clusters') ) ) {
  cat("figures/T.NK/chosen_clusters exists and is a directory.\n")
} else {
  cat("figures/T.NK/chosen_clusters does not exist - creating.\n")
  dir.create( file.path(wd, 'results', 'figures', 'T.NK', 'chosen_clusters') )
}
if (dir.exists(file.path(wd, 'results', 'figures', 'T.NK', 'UMAP') ) ) {
  cat("figures/T.NK/UMAP exists and is a directory.\n")
} else {
  cat("figures/T.NK/UMAP does not exist - creating.\n")
  dir.create( file.path(wd, 'results', 'figures', 'T.NK', 'UMAP') )
}
if (dir.exists(file.path(wd, 'results', 'figures', 'T.NK', 'Violin') ) ) {
  cat("figures/T.NK/Violin exists and is a directory.\n")
} else {
  cat("figures/T.NK/Violin does not exist - creating.\n")
  dir.create( file.path(wd, 'results', 'figures', 'T.NK', 'Violin') )
}

if (dir.exists(file.path(wd, 'results', 'calculations', 'clone_cd') ) ) {
  cat("calculations/clone_cd exists and is a directory.\n")
} else {
  cat("calculations/clone_cd does not exist - creating.\n")
  dir.create( file.path(wd, 'results', 'calculations', 'clone_cd') )
}
if (dir.exists(file.path(wd, 'results', 'figures', 'clone_cd') ) ) {
  cat("figures/clone_cd exists and is a directory.\n")
} else {
  cat("figures/clone_cd does not exist - creating.\n")
  dir.create( file.path(wd, 'results', 'figures', 'clone_cd') )
}

###################################################
# variables imported into various analysis stages #
###################################################
# project_name <- "paley_neurosarcoid_2020_08"

# project_contrasts <- list(
#   c( "case", "control")
# ) # assume each is in denominator first, numerator second order

dims_max = 30

min_cell_cutoff = 1
min_feature_cutoff = 0
max_feature_cutoff = 2500
max_mito_cutoff = 11.0

jackstraw_replicates = 500

pc_neighbor_dims = 25

mitochondrial_gene_string = "^MT-"

expected_doublet_fraction = 0.035

###############################
# color blind friendly colors #
###############################
# from colorbrewer2.org       #
###############################
colorBlindPalette <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

blues5Palette <- c( '#ffffcc', '#a1dab4', '#41b6c4', '#2c7fb8', '#253494' )
greens5Palette <- c( '#ffffcc', '#c2e699', '#78c679', '#31a354', '#006837' )
purples5Palette <- c( '#feebe2', '#fbb4b9', '#f768a1', '#c51b8a', '#7a0177' )

reds5Palette <- c( '#ffffb2', '#f3cc5c', '#fd8d3c', '#f03b20', '#bd0026' )
reds4Palette <- c( '#fef0d9', '#fdcc8a', '#fc8d59', '#d7301f' )
reds3Palette <- c( '#fee0d2', '#fc9272', '#de2d26' )

colorBlindPrintSafeDivergent1 <- c( "#a6611a", "#dfc27d", "#f5f5f5", "#80cdc1", "#018571" )
colorBlindPrintSafeDivergent2 <- c( "#7b3294", "#c2a5cf", "#f7f7f7", "#a6dba0", "#008837" )
colorBlindPrintSafeDivergent3 <- c( "#ca0020", "#f4a582", "#f7f7f7", "#92c5de", "#0571b0" )
colorBlindPrintSafeDivergent4 <- c( "#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6" )

colorBlindCopySafeSequential1 <- c( "#f0f0f0", "#bdbdbd", "#636363" )
colorBlindCopySafeSequential2 <- c( "#fef0d9", "#fdcc8a", "#fc8d59", "#d7301f" )

cluster.colors <- c(brewer.pal(12, 'Paired'), brewer.pal(8, 'Dark2'))
tissue.colors <- c( 'Eye' = 'blue' , 'Blood' = 'red' )
VUNIU.colors <- c('Viral' = 'purple', 'NIU' = 'yellow')

####################
# ggplot modifiers #
####################
gg_bigger_texts = theme(
  axis.title = element_text( size = 22 ),
  axis.text = element_text( size = 20 ),
  legend.text = element_text( size = 14 ),
  legend.title = element_text( size = 15 ),
  plot.title = element_text( size = 22 ),
  strip.text.x = element_text( size = 17, margin = margin( b = 5, t = 5 ) ),
  strip.text.y = element_text( size = 15 )
)

gg_multiplot_texts = theme(
  axis.title = element_text( size = 20 ),
  axis.text = element_text( size = 18 ),
  legend.text = element_text( size = 12 ),
  legend.title = element_text( size = 13 ),
  plot.title = element_text( size = 20 ),
  strip.text.x = element_text( size = 16, margin = margin( b = 5, t = 5 ) ),
  strip.text.y = element_text( size = 15 )
)

gg_no_legend = theme(
  legend.position='none'
)

gg_no_grid = theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

gg_no_x_grid = theme(
  panel.grid.major.x = element_blank() )

gg_no_y_grid = theme(
  panel.grid.major.y = element_blank() )

gg_center_title = theme(
  plot.title = element_text( hjust = 0.5 )
)

gg_no_x_label = theme(
  axis.title.x = element_blank()
)

gg_no_y_label = theme(
  axis.title.y = element_blank()
)

gg_angled_x_text = theme (
  axis.text.x = element_text( angle = 45, vjust = 1, hjust = 1, color = 'black' )
)

#####################
# DoMultiBarHeatmap #
#####################
suppressPackageStartupMessages({
  library(rlang)
})

DoMultiBarHeatmap <- function (object, 
                               features = NULL, 
                               cells = NULL, 
                               group.by = "ident", 
                               additional.group.by = NULL, 
                               additional.group.sort.by = NULL, 
                               cols.use = NULL,
                               group.bar = TRUE, 
                               disp.min = -2.5, 
                               disp.max = NULL, 
                               slot = "scale.data", 
                               assay = NULL, 
                               label = TRUE, 
                               size = 5.5, 
                               hjust = 0, 
                               angle = 45, 
                               raster = TRUE, 
                               draw.lines = TRUE, 
                               lines.width = NULL, 
                               group.bar.height = 0.02, 
                               combine = TRUE) 
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  
  if (!is.null(additional.group.sort.by)) {
    if (any(!additional.group.sort.by %in% additional.group.by)) {
      bad.sorts <- additional.group.sort.by[!additional.group.sort.by %in% additional.group.by]
      additional.group.sort.by <- additional.group.sort.by[additional.group.sort.by %in% additional.group.by]
      if (length(x = bad.sorts) > 0) {
        warning("The following additional sorts were omitted as they were not a subset of additional.group.by : ", 
                paste(bad.sorts, collapse = ", "))
      }
    }
  }
  
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, 
                                                             slot = slot)[features, cells, drop = FALSE])))
  
  object <- suppressMessages(expr = StashIdent(object = object, 
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by[!additional.group.by %in% group.by])]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    if (!is_null(additional.group.by)) {
      additional.group.use <- additional.group.by[additional.group.by!=i]  
      if (!is_null(additional.group.sort.by)){
        additional.sort.use = additional.group.sort.by[additional.group.sort.by != i]  
      } else {
        additional.sort.use = NULL
      }
    } else {
      additional.group.use = NULL
      additional.sort.use = NULL
    }
    
    group.use <- groups.use[, c(i, additional.group.use), drop = FALSE]
    
    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }  
    }
    
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) * 
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(rep(x = levels(x = group.use[[i]]), times = lines.width))
      group.levels <- list()
      group.levels[[i]] = levels(x = group.use[[i]])
      for (j in additional.group.use) {
        group.levels[[j]] <- levels(x = group.use[[j]])
        placeholder.groups[[j]] = NA
      }
      
      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells
      
      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells
      
      group.use <- rbind(group.use, placeholder.groups)
      
      for (j in names(group.levels)) {
        group.use[[j]] <- factor(x = group.use[[j]], levels = group.levels[[j]])
      }
      
      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }
    
    order_expr <- paste0('order(', paste(c(i, additional.sort.use), collapse=','), ')')
    group.use = with(group.use, group.use[eval(parse(text=order_expr)), , drop=F])
    
    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
                                     disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                                     cell.order = rownames(x = group.use), group.by = group.use[[i]])
    
    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))) {
        if (colname == i) {
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }
        
        # Default
        cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))))  
        
        #Overwrite if better value is provided
        if (!is_null(cols.use[[colname]])) {
          req_length = length(x = levels(group.use))
          if (length(cols.use[[colname]]) < req_length){
            warning("Cannot use provided colors for ", colname, " since there aren't enough colors.")
          } else {
            if (!is_null(names(cols.use[[colname]]))) {
              if (all(levels(group.use[[colname]]) %in% names(cols.use[[colname]]))) {
                cols[[colname]] <- as.vector(cols.use[[colname]][levels(group.use[[colname]])])
              } else {
                warning("Cannot use provided colors for ", colname, " since all levels (", paste(levels(group.use[[colname]]), collapse=","), ") are not represented.")
              }
            } else {
              cols[[colname]] <- as.vector(cols.use[[colname]])[c(1:length(x = levels(x = group.use[[colname]])))]
            }
          }
        }
        
        # Add white if there's lines
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)  
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(cols[[colname]], "#FFFFFF")
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])
        
        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)
        
        plot <- suppressMessages(plot + 
                                   annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
                                   annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                   coord_cartesian(ylim = c(0, y.max), clip = "off")) 
       
        if ((colname == i) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% pbuild$layout$panel_params[[1]]$x$break_positions()
          group.use$x <- x.divs
          label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]],
                                FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                    label.x.pos)
          plot <- plot + geom_text(stat = "identity", 
                                   data = label.x.pos, aes_string(label = "group", 
                                                                  x = "label.x.pos"), y = y.max + y.max * 
                                     0.03 * 0.5, angle = angle, hjust = hjust, 
                                   size = size)
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                                                                   y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) * 
                                                                     size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}

