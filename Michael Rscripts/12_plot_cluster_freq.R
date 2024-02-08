##############################################
# ANNOTATION AND SUBSETTING OF SEURAT OBJECT #
##############################################

#load libaries
library( ggplot2 )
library( dplyr ) 
library( tidyr )
# library( readr )
# library( purrr )
library( tibble )
library( stringr )
library( Seurat )
library( sctransform )
# library( DoubletFinder )
library(RColorBrewer)

#Increase memory usage
options(future.globals.maxSize = 60000 * 1024^2) #60GB

# Set Working Directory
wd <- getwd()
# source
source( file = file.path(wd, "src", "shared_r_functions.R" ) )

# Load the Seurat Object 'scd' ("Single Cell Data")
rds.path <- file.path( object_path , "scd.annotated.rds" )
scd <- readRDS(rds.path)
scd # Sanity check
head( Idents(object = scd) ) # Sanity check

# Set variables
section <- 'infection.by.celltype'
n <- length( unique( scd$lin_1 ) )

# Table frequencies
counts <- table( scd$lin_1, scd$Infection )
freq <- round(counts/rowSums(counts), digits = 4) #normalize to percentages within each cell state

# Set color scheme base for plots
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Save calulations
write.csv( as.matrix(counts), file = file.path(wd, 'results', 'calculations', paste0(section, '.counts.per.sample.matrix.csv') ), quote = FALSE, row.names = TRUE )
write.csv( as.matrix(freq), file = file.path(wd, 'results', 'calculations', paste0(section,'.freq.per.sample.matrix.csv') ), quote = FALSE, row.names = TRUE )
# Convert Cluster Frequencies Table to Dataframe for plotting
counts <- as.data.frame(counts)
colnames(counts) <- c('lin_1', 'Infection', 'Counts') 
freq <- as.data.frame(freq)
colnames(freq) <- c('lin_1', 'Infection', 'Freq') 
# Save calulations
write.csv( as.data.frame(counts), file = file.path(wd, 'results', 'calculations', paste0(section,'.counts.per.sample.df.csv') ), quote = FALSE, row.names = TRUE )
write.csv( as.data.frame(freq), file = file.path(wd, 'results', 'calculations', paste0(section,'.freq.per.sample.df.csv') ), quote = FALSE, row.names = TRUE )

# Print Cluster Frequencies in Stacked Bar Plot
file.name <- paste(section,'freq','per','sample','pdf', sep='.')
pdf( file = file.path(wd, 'results', 'calculations', file.name ) )
ggplot(counts, aes( fill = lin_1, y = Counts, x = Infection ) ) + # Samples are ordered by njoin
  geom_bar(position="fill", stat="identity") +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual( values = sample(col_vector, n), aesthetics = "fill" ) + # Using 'n' defined above
  theme(panel.background=element_rect(fill="white"), # background = white
        axis.text.x = element_text(angle=45, hjust = 1,vjust=1,size = 12,face = "bold"),
        axis.text.y = element_text(size = 12,face = "bold"))
dev.off()
cat( paste0(section, " ", "section complete.\n") )
rm( section, n, counts, freq )

gc()


