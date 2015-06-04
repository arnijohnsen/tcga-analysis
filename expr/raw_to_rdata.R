library(data.table)

# Define cancer type, raw and parsed data directories --------------------------
cancer_type <- "brca"
raw_data_dir    <- "/share/scratch/arj32/raw_data/"
parsed_data_dir <- "/share/scratch/arj32/parsed_data/"
source_file_path<- "/expr/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/"
source_file_dir <- paste(raw_data_dir, cancer_type, source_file_path, sep="")
output_dir      <- paste(parsed_data_dir, cancer_type, "/expr/",  sep="")

# Load file and sample information ---------------------------------------------
file_sample_map <- fread(paste(raw_data_dir, cancer_type,
                               "/expr/FILE_SAMPLE_MAP.txt", sep=""))
setnames(file_sample_map, c("filename", "barcode"))
# use only genes_normalized_results
file_sample_map <- file_sample_map[grep("genes_normalized_results", filename)]

sample_list <- fread(paste(parsed_data_dir, cancer_type,
                               "/info/expr_participants.txt", sep=""))

n <- dim(sample_list)[1]

# Read one file with all information, use to create data.tables ----------------
cat("Initiating tables..\n")
idx_cancer <- which(!is_na(sample_list$cancer_barcode))[1]
idx_normal <- which(!is_na(sample_list$normal_barcode))[1]
cancer_barcode <- sample_list$cancer_barcode[idx_cancer]
normal_barcode <- sample_list$normal_barcode[idx_normal]
cancer_file <- file_sample_map[barcode==cancer_barcode, filename]
normal_file <- file_sample_map[barcode==normal_barcode, filename]

cancer_wide <- fread(paste(source_file_dir, cancer_file, sep=""))
normal_wide <- fread(paste(source_file_dir, normal_file, sep=""))

cancer_wide[,normalized_count:=NULL]
normal_wide[,normalized_count:=NULL]

setnames(cancer_wide, c("gene"))
setnames(normal_wide, c("gene"))

cancer_wide[,gene:=gsub("\\|.*", "",gene)]
normal_wide[,gene:=gsub("\\|.*", "",gene)]

# Read all files, drop everything except beta and append as a new column--------
for(i in 1:n){
  cat("Reading files for participant", i, "of", n, "\n")
  cancer_barcode <- sample_list$cancer_barcode[i]
  normal_barcode <- sample_list$normal_barcode[i]
  cancer_file <- file_sample_map[barcode==cancer_barcode, filename]
  normal_file <- file_sample_map[barcode==normal_barcode, filename]
  if(!is_na(cancer_barcode)){
    tmp <- fread(paste(source_file_dir, cancer_file, sep=""), select=2)
    cancer_wide[,c(cancer_barcode):=tmp]
  }
  if(!is.na(normal_barcode)){
    tmp <- fread(paste(source_file_dir, normal_file, sep=""), select=2)
    normal_wide[,c(normal_barcode):=tmp]
  }
}

# Assign systematic names to data frames and save ------------------------------
assign(paste(cancer_type, "_expr_cancer", sep=""), cancer_wide)
assign(paste(cancer_type, "_expr_normal", sep=""), normal_wide)
save(list = paste(cancer_type, "_expr_cancer", sep=""),
     file = paste(output_dir, cancer_type, "_expr_cancer.Rdata", sep=""))
save(list = paste(cancer_type, "_expr_normal", sep=""),
     file = paste(output_dir, cancer_type, "_expr_normal.Rdata", sep=""))
