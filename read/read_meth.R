library(data.table)
library(yaml)    # Read config
library(utils)   # Progress bar
library(R.utils) # Verbose

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

# Define absolute directories --------------------------------------------------
cat(verbose, "Defining directories to read from")
file_filter <- ".*"
src_file_path <- switch(
  cancer_type,
  ov = "/DNA_Methylation/JHU_USC__HumanMethylation27/Level_3/",
  "/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"
)
src_file_dir <- paste(
  raw_data_dir, cancer_type, "/", data_type, src_file_path, sep = ""
)
output_dir <- paste(parsed_data_dir, cancer_type, "/", data_type, "/",  sep="")

# Generate list of files to read and barcodes of used samples ------------------
cat(verbose, "Generating list of files to read")
file_sample_map <- fread(
  paste(
    raw_data_dir, cancer_type, "/", data_type, "/FILE_SAMPLE_MAP.txt", sep = ""
  )
)
setnames(file_sample_map, c("filename", "barcode"))
file_sample_map <- file_sample_map[grep(file_filter, filename)]
setkey(file_sample_map, barcode)
all_participants <- fread(
  paste(
    parsed_data_dir, cancer_type, "/info/participants.txt", sep = ""
  )
)
participants <- unlist(
  all_participants[, grep(data_type, names(all_participants)), with = F]
)
participants <- participants[!is.na(participants)]
exit(verbose)

# Create first columns of tables, which contain names of genes, probes, etc. ---
enter(verbose, "Reading source files")
cat(verbose, "Reading from ", src_file_dir)
cat(verbose, "Initiating tables")
cancer_table <- fread(
  paste(
    src_file_dir, file_sample_map[participants[1]]$filename, sep = ""
  ), 
  skip = 1, showProgress = FALSE
)
cancer_table[,Beta_value:=NULL]
setnames(cancer_table, c("probe", "gene", "chrom", "coord"))
normal_table <- copy(cancer_table)

# Read files one by one and add to table ---------------------------------------
cat(verbose, "Progress:")
n <- length(participants)
pb <- txtProgressBar(max = n, style = 3)
for(i in 1:n){
  tmp <- fread(
    paste(
      src_file_dir, file_sample_map[participants[i]]$filename, sep = ""
    ),
    skip = 1, select = 2, showProgress = FALSE
  )
  if (substr(participants[i], 14, 15) == "01") {
    cancer_table[,c(participants[i]):=tmp]
  } else {
    normal_table[,c(participants[i]):=tmp]
  }
  setTxtProgressBar(pb, i)
}
newline(verbose)
exit(verbose)

# Save to .RDS file ------------------------------------------------------------
enter(verbose, "Saving files")
cat(
  verbose, "Saving to cancer file to ",
  paste(output_dir, cancer_type, "_", data_type, "_cancer.Rds", sep = "")
)
saveRDS(
  cancer_table,
  file = paste(
    output_dir, cancer_type, "_", data_type, "_cancer_full.Rds", sep = ""
  )
)
cat(
  verbose, "Saving to normal file to ",
  paste(output_dir, cancer_type, "_", data_type, "_normal.Rds", sep = "")
)
saveRDS(
  normal_table,
  file = paste(
    output_dir, cancer_type, "_", data_type, "_normal_full.Rds", sep = ""
  )
)
exit(verbose)
