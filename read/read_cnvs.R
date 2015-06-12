library(data.table)
library(yaml)    # Read config
library(utils)   # Progress bar
library(R.utils) # Verbose
library(CNTools) # Convert from long to wide

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
data_type       <- "cnvs"

# Define absolute directories --------------------------------------------------
cat(verbose, "Defining directories to read from")
file_filter  <- "nocnv_hg19"
src_file_path <- "/CNV_SNP_Array/BI__Genome_Wide_SNP_6/Level_3/"
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

# Create empty lists where data frames are added -------------------------------
enter(verbose, "Reading source files")
cat(verbose, "Creating empty lists")
cnvs_list <- list()

# Reading files one by one and adding to list ----------------------------------
cat(verbose, "Reading files from", src_file_dir)
cat(verbose, "Progress:")
n <- length(participants)
pb <- txtProgressBar(max = n, style = 3)
for(i in 1:n){
  tmp <- fread(
    paste(
      src_file_dir, file_sample_map[participants[i]]$filename, sep = ""
    )
  )
  tmp$Sample <- participants[i]
  cnvs_list[[participants[i]]] <- tmp
  setTxtProgressBar(pb, i)
}
newline(verbose)
exit(verbose)

# Binding data to one data frame and transforming with CNTools -----------------
enter(verbose, "Converting data to wide format")
cat(verbose, "Binding list to one data.table")
cnvs_bound <- rbindlist(cnvs_list)
setnames(
  cnvs_bound,
  c("ID", "chrom", "loc.start", "loc.end", "n.probes", "seg.mean")
)

cat(verbose, "Converting from long to wide format")
gene_info <- fread(
  paste(
    parsed_data_dir, cancer_type, "/info/hg19geneinfo.txt", sep=""
  )
)
cnvs_seg <- CNSeg(cnvs_bound)
cnvs_rs  <- getRS(
  cnvs_seg,
  by = "gene", imput = F, XY = T, what = "mean", geneMap = gene_info
)
cnvs_cancer <-data.table(rs(cnvs_rs))
setnames(cnvs_cancer, 4, "gene")
exit(verbose)

# Save to .RDS file ------------------------------------------------------------
enter(verbose, "Saving files")
cat(
  verbose, "List format to ",
  paste(output_dir, cancer_type, "_cnvl_cancer.Rds", sep = "")
)
saveRDS(
  cnvs_list,
  file = paste(output_dir, cancer_type, "_cnvl_cancer.Rds", sep = "")
)
cat(
  verbose, "Wide format to ",
  paste(output_dir, cancer_type, "_cnvs_cancer.Rds", sep = "")
)
saveRDS(
  cnvs_cancer,
  file = paste(output_dir, cancer_type, "_cnvs_cancer.Rds", sep = "")
)
exit(verbose)
