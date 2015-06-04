library(data.table)

# Define cancer type, raw and parsed data directories --------------------------
cancer_type <- "brca"
raw_data_dir    <- "/share/scratch/arj32/raw_data/"
parsed_data_dir <- "/share/scratch/arj32/parsed_data/"
source_file_path<- "/meth/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"
source_file_dir <- paste(raw_data_dir, cancer_type, source_file_path, sep="")
output_dir      <- paste(parsed_data_dir, cancer_type, "/meth/",  sep="")

# Load file and sample information ---------------------------------------------
file_sample_map <- fread(paste(raw_data_dir, cancer_type,
                               "/meth/FILE_SAMPLE_MAP.txt", sep=""))
setnames(file_sample_map, c("filename", "barcode"))
sample_list <- fread(paste(parsed_data_dir, cancer_type,
                               "/info/meth_participants.txt", sep=""))

n <- dim(sample.list)[1]

# Read one file with all information, use to create data_tables ----------------
cat("Initiating tables..\n")
idx_cancer <- which(!is_na(sample_list$cancer_barcode))[1]
idx_normal <- which(!is_na(sample_list$normal_barcode))[1]
cancer_barcode <- sample_list$cancer_barcode[idx_cancer]
normal_barcode <- sample_list$normal_barcode[idx_normal]
cancer_file <- file_sample_map[barcode==cancer_barcode, filename]
normal_file <- file_sample_map[barcode==normal_barcode, filename]

cancer_wide <- fread(paste(source_file_dir, cancer_file, sep=""), skip=1)
normal_wide <- fread(paste(source_file_dir, normal_file, sep=""), skip=1)

cancer_wide[,Beta_value:=NULL]
normal_wide[,Beta_value:=NULL]

setnames(cancer_wide, c("probe", "gene", "chrom", "coord"))
setnames(normal_wide, c("probe", "gene", "chrom", "coord"))

# Read all files, drop everything except beta and append as a new column--------
for(i in 1:n){
  cat("Reading files for participant", i, "of", n, "\n")
  cancer_barcode <- sample_list$cancer_barcode[i]
  normal_barcode <- sample_list$normal_barcode[i]
  cancer_file <- file_sample_map[barcode==cancer_barcode, filename]
  normal_file <- file_sample_map[barcode==normal_barcode, filename]
  if(!is_na(cancer_barcode)){
    tmp <- fread(paste(source_file_dir, cancer_file, sep=""), skip=1, select=2)
    cancer_wide[,c(cancer_barcode):=tmp]
  }
  if(!is.na(normal_barcode)){
    tmp <- fread(paste(source_file_dir, normal_file, sep=""), skip=1, select=2)
    normal_wide[,c(normal_barcode):=tmp]
  }
}

# Assign systematic names to data frames and save ------------------------------
assign(paste(cancer_type, "_meth_cancer", sep=""), cancer_wide)
assign(paste(cancer_type, "_meth_normal", sep=""), normal_wide)
save(list = paste(cancer_type, "_meth_cancer", sep=""),
     file = paste(output_dir, cancer_type, "_meth_cancer.Rdata", sep=""))
save(list = paste(cancer_type, "_meth_normal", sep=""),
     file = paste(output_dir, cancer_type, "_meth_normal.Rdata", sep=""))
