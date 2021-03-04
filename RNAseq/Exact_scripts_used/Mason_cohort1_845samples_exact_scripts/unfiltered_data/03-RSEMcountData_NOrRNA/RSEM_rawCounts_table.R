library(tximport)
library(tidyverse)

work_dir="/path/to/COVIRT_all_data_rRNA_removed/processing_scripts/03-RSEMcountData"
counts_dir="/path/to/COVIRT_all_data_rRNA_removed/03-RSEMcountData"

setwd(file.path(work_dir))

### Pull in sample metadata ###
study <- read.csv(Sys.glob(file.path(work_dir,"samples.txt")), header = FALSE, row.names = 1, stringsAsFactors = TRUE)


##### Import Data
files <- list.files(file.path(counts_dir),pattern = ".genes.results", full.names = TRUE)
# reorder the genes.results files to match the ordering of the samples in the metadata file
files <- files[sapply(rownames(study), function(x)grep(x, files, value=FALSE, fixed=TRUE))]
names(files) <- rownames(study)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)


##### Export unnormalized counts table
setwd(file.path(counts_dir))
write.csv(txi.rsem$counts,file='RSEM_Unnormalized_Counts.csv')


## print session info ##
print("Session Info below: ")
print("")
sessionInfo()
