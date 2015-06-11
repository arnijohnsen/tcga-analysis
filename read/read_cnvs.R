library(data.table)

# Define cancer type, raw and parsed data directories --------------------------
cancer_type <- "brca"
raw_data_dir    <- "/share/scratch/arj32/raw_data/"
parsed_data_dir <- "/share/scratch/arj32/parsed_data/"
source_file_path<- "/cnvs/CNV_SNP_Array/BI__Genome_Wide_SNP_6/Level_3/"
source_file_dir <- paste(raw_data_dir, cancer_type, source_file_path, sep="")
output_dir      <- paste(parsed_data_dir, cancer_type, "/cnvs/",  sep="")

# Load file and sample information ---------------------------------------------
file_sample_map <- fread(paste(raw_data_dir, cancer_type,
                               "/cnvs/FILE_SAMPLE_MAP.txt", sep=""))
setnames(file_sample_map, c("filename", "barcode"))
# Use only nocsv_hg19
file_sample_map <- file_sample_map[grep("nocnv_hg19", filename)]

sample_list <- fread(paste(parsed_data_dir, cancer_type,
                           "/info/cnvs_participants.txt", sep=""))

n <- dim(sample_list)[1]

# Create empty lists -----------------------------------------------------------
cancer_list <- list()
normal_list <- list()

# Read files and create a list of data.tables ----------------------------------
for(i in 1:n){
  cat("Reading files for participant", i, "of", n, "\n")
  cancer_barcode <- sample_list$cancer_barcode[i]
  normal_barcode <- sample_list$normal_barcode[i]
  cancer_file <- file_sample_map[barcode==cancer_barcode, filename]
  normal_file <- file_sample_map[barcode==normal_barcode, filename]
  if(!is.na(cancer_barcode)){
    tmp <- fread(paste(source_file_dir, cancer_file, sep=""))
    tmp$Sample <- cancer_barcode
    cancer_list[[ cancer_barcode ]] <- tmp
  }
  if(!is.na(normal_barcode)){
    tmp <- fread(paste(source_file_dir, normal_file, sep=""))
    tmp$Sample <- normal_barcode
    normal_list[[ normal_barcode]] <- tmp
  }
}

# Assign systematic names to data frames and save ------------------------------
assign(paste(cancer_type, "_cnvl_cancer", sep=""), cancer_list)
assign(paste(cancer_type, "_cnvl_normal", sep=""), normal_list)
save(list = paste(cancer_type, "_cnvl_cancer", sep=""),
     file = paste(output_dir, cancer_type, "_cnvl_cancer.Rdata", sep=""))
save(list = paste(cancer_type, "_cnvl_normal", sep=""),
     file = paste(output_dir, cancer_type, "_cnvl_normal.Rdata", sep=""))
