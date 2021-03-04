print("Make STAR counts table")
print("")

work_dir="/path/to/COVIRT_all_data_rRNA_removed/processing_scripts/02-AlignedData"
align_dir="/path/to/COVIRT_all_data_rRNA_removed/02-AlignedData/STAR_counts"

setwd(file.path(work_dir))

### Pull in sample names ###
study <- read.csv(Sys.glob(file.path(work_dir,"samples.txt")), header = FALSE, row.names = 1, stringsAsFactors = TRUE)

##### Import Data
ff <- list.files(file.path(align_dir),pattern = "ReadsPerGene.out.tab", full.names = TRUE)
# Remove the first 4 lines
counts.files <- lapply( ff, read.table, skip = 4 )
# Get counts aligned to the second, reverse, strand
counts <- as.data.frame( sapply( counts.files, function(x) x[ , 4 ] ) )
# Add column and row names
colnames(counts) <- rownames(study)
row.names(counts) <- counts.files[[1]]$V1


##### Export unnormalized counts table
setwd(file.path(align_dir))
write.csv(counts,file='STAR_Unnormalized_Counts.csv')


## print session info ##
print("Session Info below: ")
print("")
sessionInfo()
