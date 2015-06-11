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
data_type       <- "muta"

# Define absolute directories --------------------------------------------------
cat(verbose, "Defining directories to read from")
file_filter <- ".*"
src_file_path <- paste(
  "/Somatic_Mutations/",
  switch(cancer_type,
    brca = "WUSM__IlluminaGA",
    ov   = "BCM__SOLiD",
    paad = "BI__IlluminaGA",
    prad = "BI__IlluminaGA"
  ),
  "_DNASeq_curated/Level_2/", sep = "")
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

# Read from only data file  ----------------------------------------------------
exit(verbose)
enter(verbose, "Reading source files")
cat(verbose, "Reading from ", src_file_dir)

mutation <- fread(paste(src_file_dir, file_sample_map[1,filename], sep=""),
                  select=c(1,5,6,7,9,16,17))
setnames(mutation, c("gene", "chrom", "start", "end",
                     "type", "cancer_barcode", "normal_barcode"))
exit(verbose)

# Save to .RDS file ------------------------------------------------------------
enter(verbose, "Saving files")
cat(
  verbose, "Saving to cancer file to ",
  paste(output_dir, cancer_type, "_", data_type, "_cancer.Rds", sep = "")
)
saveRDS(
  mutation,
  file = paste(
    output_dir, cancer_type, "_", data_type, "_cancer.Rds", sep = ""
  )
)
exit(verbose)
