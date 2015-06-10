library(data.table)
library(yaml)

config <- yaml.load_file("/share/scratch/arj32/tcga-analysis/config.yml")

raw_data_dir    <- config$data_dirs$raw_data
parsed_data_dir <- config$data_dirs$parsed_data
cancer_type     <- config$cancer_type

barcode_regex <- paste("^TCGA",         # TCGA Project
                       "[A-Z0-9]{2}",   # Tissue source site
                       "[A-Z0-9]{4}",   # Participant
                       "[01]{2}[A-Z]",  # Sample and vial
                       "[0-9]{2}[A-Z]", # Porion and analyte
                       "[A-Z0-9]{4}",   # Plate
                       "[A-Z0-9]{2}$",  # Center
                       sep = "-")
data_types <- c("cnvs", "expr", "meth", "mirn", "muta")
file_filter <- c("nocnv_hg19",                 # cnvs
                 "genes\\.normalized_results", # expr
                 ".*",                         # meth
                 "mirna\\.quantification",     # mirn
                 ".*")                         # muta
names(file_filter) <- data_types

barcode_list <- list()

for(data_type in data_types){
  if (data_type == "muta"){
    # Specific method for muta -------------------------------------------------
    filename <- fread(paste(raw_data_dir, cancer_type, "/", data_type,
                            "/FILE_SAMPLE_MAP.txt", sep=""))[1,filename]
    file_sample_map <- fread(paste(raw_data_dir, cancer_type, "/", data_type,
                                   "/Somatic_Mutations/",
                                   "WUSM__IlluminaGA_DNASeq_curated/Level_2/",
                                   filename, sep=""), select=c(16))
    file_sample_map <- unique(file_sample_map)
    file_sample_map[, type:="muta_c"]
    barcode_list[[data_type]] <- file_sample_map
  } else {
    # General method for cnvs, expr, meth and mirn -----------------------------

    # Read FILE_SAMPLE_MAP.txt and set names -----------------------------------
    file_sample_map <- fread(paste(raw_data_dir, cancer_type, "/", data_type,
                                   "/FILE_SAMPLE_MAP.txt", sep = ""))

    setnames(file_sample_map, c("filename", "barcode"))
    file_sample_map <- file_sample_map[grep(file_filter[data_type], filename)]
    file_sample_map[, filename:=NULL]
    file_sample_map[, type:=paste(data_type,
                                  ifelse(substr(barcode, 14,14) == "0",
                                         "c", "n"), sep = "_")]

    barcode_list[[data_type]] <- file_sample_map

  }
}

barcodes_bound <- rbindlist(barcode_list)
barcodes_bound <- barcodes_bound[type!="cnvs_n"]

# Remove all samples which don't match a barcode regex -----------------------
barcodes_bound <- barcodes_bound[grep(barcode_regex, barcode)]
barcodes_bound[, participant := substr(barcode, 1, 12)]
setkey(barcodes_bound, participant)

# Define function to aggregate and then aggregate ----------------------------
select_barcode <- function(x){
  if (length(x) > 1){
      warning(paste("Unhandled multiples:",
                    paste(x, collapse=" ")))
  }
  return(x[1])
}

casted_barcodes <- dcast.data.table(barcodes_bound, participant ~ type,
                                    value.var = "barcode",
                                    fun.aggregate = select_barcode)
