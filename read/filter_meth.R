library(data.table)
library(yaml)    # Read config
library(utils)   # Progress bar
library(R.utils) # Verbose
library(stringr) # String parsing

# ONLY RUN ON MACHINES WITH PLENTY OF RAM (4 GB IS NOT ENOUGH) -----------------

# Define verbose and verbose level
verbose <- Verbose(threshold = -1)

# Read base directories from config.yml, define cancer and data type -----------
ruler(verbose)
enter(verbose, "Preparation")
cat(verbose, "Reading config.yml")
config <- yaml.load_file("/share/scratch/arj32/tcga-analysis/config.yml")
raw_data_dir    <- config$data_dirs$raw_data
parsed_data_dir <- config$data_dirs$parsed_data
cancer_type     <- config$cancer_type
data_type       <- "meth"
output_dir <- paste(parsed_data_dir, cancer_type, "/", data_type, "/",  sep="")
exit(verbose)

# Define function to count how many intervals a certain probe appreas in -------
count_intervals <- function(x){
  i02 <- x < 0.2
  i04 <- x < 0.4
  i06 <- x < 0.6
  i08 <- x < 0.8
  intervals <- c(
    sum(i02, na.rm = T),
    sum(!i02 & i04, na.rm = T),
    sum(!i04 & i06, na.rm = T),
    sum(!i06 & i08, na.rm = T),
    sum(!i08, na.rm = T)
  )
  return(sum(intervals != 0))
}

enter(verbose, "Reading cancer file")
meth_cancer_full <- readRDS(
  paste(
    output_dir, cancer_type, "_", data_type, "_cancer_full.Rds", sep = ""
  )
)
exit(verbose)
enter(verbose, "Finding good probes")
cat(verbose, "Setting key")
setkey(meth_cancer_full, probe)
cat(verbose, "Applying function")
count <- apply(meth_cancer_full[,5:728, with = F], 1, count_intervals)
good_probes <- meth_cancer_full[count > 2]$probe
exit(verbose)

enter(verbose, "Saving cancer file")
cat(
  verbose, "Saving to cancer file to ",
  paste(output_dir, cancer_type, "_", data_type, "_cancer.Rds", sep = "")
)
saveRDS(
  meth_cancer_full[good_probes, 1:728, with=F],
  file = paste(
    output_dir, cancer_type, "_", data_type, "_cancer.Rds", sep = ""
  )
)
cat(verbose, "Deleting cancer file")
rm(meth_cancer_full)
gc()
exit(verbose)

enter(verbose, "Reading normal file")
meth_normal_full <- readRDS(
  paste(
    output_dir, cancer_type, "_", data_type, "_normal_full.Rds", sep = ""
  )
)
setkey(meth_normal_full, probe)
exit(verbose)
enter(verbose, "Saving normal file")
cat(
  verbose, "Saving to normal file to ",
  paste(output_dir, cancer_type, "_", data_type, "_normal.Rds", sep = "")
)
saveRDS(
  meth_normal_full[good_probes, 1:100, with=F],
  file = paste(
    output_dir, cancer_type, "_", data_type, "_normal.Rds", sep = ""
  )
)
exit(verbose)
