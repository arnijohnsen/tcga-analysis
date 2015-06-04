library(data.table)

# Define cancer type, raw and parsed data directories --------------------------
cancer_type <- "brca"
raw_data_dir    <- "/share/scratch/arj32/raw_data/"
parsed_data_dir <- "/share/scratch/arj32/parsed_data/"
source_file_path<- "/mirn/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3/"
source_file_dir <- paste(raw_data_dir, cancer_type, source_file_path, sep="")
output_dir      <- paste(parsed_data_dir, cancer_type, "/mirn/",  sep="")

# Load file and sample information ---------------------------------------------
file_sample_map <- fread(paste(raw_data_dir, cancer_type,
                               "/mirn/FILE_SAMPLE_MAP.txt", sep=""))
setnames(file_sample_map, c("filename", "barcode"))
# use only genes_normalized_results
file_sample_map <- file_sample_map[grep("mirna.quantification", filename)]

sample_list <- fread(paste(parsed_data_dir, cancer_type,
                               "/info/mirn_participants_txt", sep=""))

sample_list[normal_barcode == ""]$normal_barcode <- NA

n <- dim(sample_list)[1]

# Read one file with all information, use to create data.tables ----------------
cat("Initiating tables..\n")
idx_cancer <- which(!is_na(sample_list$cancer_barcode))[1]
idx_normal <- which(!is_na(sample_list$normal_barcode))[1]
cancer_barcode <- sample_list$cancer_barcode[idx_cancer]
normal_barcode <- sample_list$normal_barcode[idx_normal]
cancer_file <- file_sample_map[barcode==cancer_barcode, filename]
normal_file <- file_sample_map[barcode==normal_barcode, filename]

cat(paste(source_file_dir, cancer_file, sep=""), "\n")
cat(paste(source_file_dir, normal_file, sep=""), "\n")

cancer_wide <- fread(paste(source_file_dir, cancer_file, sep=""),select=1)
normal_wide <- fread(paste(source_file_dir, normal_file, sep=""),select=1)

setnames(cancer_wide, c("mirna"))
setnames(normal_wide, c("mirna"))

# Read all files, drop everything except beta and append as a new column--------
for(i in 1:n){
  cat("Reading files for participant", i, "of", n, "\n")
  cancer_barcode <- sample_list$cancer_barcode[i]
  normal_barcode <- sample_list$normal_barcode[i]
  cancer_file <- file_sample_map[barcode==cancer_barcode, filename]
  normal_file <- file_sample_map[barcode==normal_barcode, filename]
  if(!is.na(cancer_barcode)){
    tmp <- fread(paste(source_file_dir, cancer_file, sep=""), select=3)
    cancer_wide[,c(cancer_barcode):=tmp]
  }
  if(!is.na(normal_barcode)){
    tmp <- fread(paste(source_file_dir, normal_file, sep=""), select=3)
    normal_wide[,c(normal_barcode):=tmp]
  }
}

# Assign systematic names to data frames and save ------------------------------
assign(paste(cancer_type, "_mirn_cancer", sep=""), cancer_wide)
assign(paste(cancer_type, "_mirn_normal", sep=""), normal_wide)
save(list = paste(cancer_type, "_mirn_cancer", sep=""),
     file = paste(output_dir, cancer_type, "_mirn_cancer.Rdata", sep=""))
save(list = paste(cancer_type, "_mirn_normal", sep=""),
     file = paste(output_dir, cancer_type, "_mirn_normal.Rdata", sep=""))
