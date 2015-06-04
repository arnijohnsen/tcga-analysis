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
database_info <- fread(paste(source_file_dir, file_sample_map[1,filename], sep=""),
                     select=c(16,17))
setnames(database_info, c("cancer_barcode", "normal_barcode"))
database_info[,participant:=substr(cancer_barcode,1,12)]
setcolorder(database_info, c("participant", "cancer_barcode", "normal_barcode"))

# Filter unique samples
setkey(database_info, participant, cancer_barcode, normal_barcode)
database_info_unq <- unique(database_info)

# Assign systematic names to data frames and save ------------------------------
write.table(database_info_unq, paste(parsed_data_dir, cancer_type,
                                     "/info/muta_participants.txt", sep=""),
            quote=F, row.names=F)
