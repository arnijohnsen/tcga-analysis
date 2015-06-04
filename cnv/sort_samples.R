library(data.table)

cancer_type <- "brca"
raw_data_dir <- "/share/scratch/arj32/raw_data/"
parsed_data_dir <- "/share/scratch/arj32/parsed_data/"

file_sample_map <- fread(paste(raw_data_dir, cancer_type, "/cnv/FILE_SAMPLE_MAP.txt", sep=""))

# use only nocnv_hg19
file_sample_map <- file_sample_map[grep("nocnv_hg19", filename)]

# dont use recurrent cancers or files with multiple barcodes
bad_barcodes <- grepl(",", file_sample_map$barcode) | (substr(file_sample_map$barcode, 14, 15) == "02")
file_sample_map <- file_sample_map[!bad_barcodes]

database_info <- data.table(participant = unique(substr(file_sample_map$barcode, 1, 12)))

# Find cancer barcode
cancer_barcode <- sapply(database_info$participant, function(x){tmp <- grep(paste(x,"-01",sep=""), file_sample_map$barcode)
                                                             if(length(tmp) > 1){
                                                               warning("Sample has more than 1 cancer")
                                                               tmp <- tmp[1]
                                                             }
                                                             if(length(tmp) > 0){
                                                               file_sample_map$barcode[tmp]
                                                             } else {
                                                               NA
                                                             }
                                                             })
database_info$cancer_barcode <- unlist(cancer_barcode)

# Find normal barcode
normal_barcode <- sapply(database_info$participant, function(x){tmp <- grep(paste(x,"-1",sep=""), file_sample_map$barcode)
                                                             if(length(tmp) > 1){
                                                               warning("Sample has more than 1 normal")
                                                               tmp <- tmp[1]
                                                             }
                                                             if(length(tmp) > 0){
                                                               file_sample_map$barcode[tmp]
                                                             } else {
                                                               NA
                                                             }
                                                             })
database_info$normal_barcode <- unlist(normal_barcode)

write.table(database_info, paste(parsed_data_dir, cancer_type, "/info/cnv_participants.txt", sep=""), quote=F, row.names=F)
