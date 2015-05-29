library(data.table)

cancer.type <- "brca"
raw.data.dir    <- "/share/scratch/arj32/raw-data/"
parsed.data.dir <- "/share/scratch/arj32/parsed-data/"

file.sample.map <- fread(paste(raw.data.dir, cancer.type, "/expr/FILE_SAMPLE_MAP.txt", sep=""))

# use only genes.normalized_results
file.sample.map <- file.sample.map[grep("genes.normalized_results", filename)]

# dont use recurrent cancers or files with multiple barcodes
bad.barcodes <- grepl(",", file.sample.map$barcode) | (substr(file.sample.map$barcode, 14, 15) == "06")
#bad.barcodes <- grepl(",", file.sample.map$barcode)
file.sample.map <- file.sample.map[!bad.barcodes]

database.info <- data.table(participant = unique(substr(file.sample.map$barcode, 1, 12)))

# Find cancer barcode
cancer.barcode <- sapply(database.info$participant, function(x){tmp <- grep(paste(x,"-01",sep=""), file.sample.map$barcode)
                                                             if(length(tmp) > 1){
                                                               warning("Sample has more than 1 cancer")
                                                               tmp <- tmp[1]
                                                             }
                                                             if(length(tmp) > 0){
                                                               file.sample.map$barcode[tmp]
                                                             } else {
                                                               NA
                                                             }
                                                             })
database.info$cancer.barcode <- unlist(cancer.barcode)

# Find normal barcode
normal.barcode <- sapply(database.info$participant, function(x){tmp <- grep(paste(x,"-1",sep=""), file.sample.map$barcode)
                                                             if(length(tmp) > 1){
                                                               warning("Sample has more than 1 normal")
                                                               tmp <- tmp[1]
                                                             }
                                                             if(length(tmp) > 0){
                                                               file.sample.map$barcode[tmp]
                                                             } else {
                                                               NA
                                                             }
                                                             })
database.info$normal.barcode <- unlist(normal.barcode)

write.table(database.info, paste(parsed.data.dir, cancer.type, "/info/expr-participants.txt", sep=""), quote=F, row.names=F)

