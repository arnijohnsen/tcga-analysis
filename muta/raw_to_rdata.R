library(data.table)

# Define cancer type, raw and parsed data directories --------------------------
cancer_type <- "brca"
raw_data_dir    <- "/share/scratch/arj32/raw_data/"
parsed_data_dir <- "/share/scratch/arj32/parsed_data/"
source_file_path<- "/muta/Somatic_Mutations/WUSM__IlluminaGA_DNASeq_curated/Level_2/"
source_file_dir <- paste(raw_data_dir, cancer_type, source_file_path, sep="")
output_dir      <- paste(parsed_data_dir, cancer_type, "/muta/",  sep="")

# Load file and sample information ---------------------------------------------
file_sample_map <- fread(paste(raw_data_dir, cancer_type,
                               "/muta/FILE_SAMPLE_MAP.txt", sep=""))
setnames(file_sample_map, c("filename", "barcode"))

# Read only file ---------------------------------------------------------------
mutation <- fread(paste(source_file_dir, file_sample_map[1,filename], sep=""),
                  select=c(1,5,6,7,9,16,17))
setnames(mutation, c("gene", "chrom", "start", "end",
                     "type", "cancer_barcode", "normal_barcode"))

# Assign systematic names to data frames and save ------------------------------
assign(paste(cancer_type, "_muta_cancer", sep=""), mutation)
save(list = paste(cancer_type, "_muta_cancer", sep=""),
     file = paste(output_dir, cancer_type, "_muta_cancer.Rdata", sep=""))
