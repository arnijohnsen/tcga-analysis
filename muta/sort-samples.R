library(data.table)

# Define cancer type, raw and parsed data directories --------------------------
cancer.type <- "brca"
raw.data.dir    <- "/share/scratch/arj32/raw-data/"
parsed.data.dir <- "/share/scratch/arj32/parsed-data/"
source.file.path<- "/muta/Somatic_Mutations/WUSM__IlluminaGA_DNASeq_curated/Level_2/"
source.file.dir <- paste(raw.data.dir, cancer.type, source.file.path, sep="")
output.dir      <- paste(parsed.data.dir, cancer.type, "/muta/",  sep="")

# Load file and sample information ---------------------------------------------
file.sample.map <- fread(paste(raw.data.dir, cancer.type,
                               "/muta/FILE_SAMPLE_MAP.txt", sep=""))
setnames(file.sample.map, c("filename", "barcode"))

# Read only file ---------------------------------------------------------------
database.info <- fread(paste(source.file.dir, file.sample.map[1,filename], sep=""), 
                     select=c(16,17))
setnames(database.info, c("cancer.barcode", "normal.barcode"))
database.info[,participant:=substr(cancer.barcode,1,12)]
setcolorder(database.info, c("participant", "cancer.barcode", "normal.barcode"))

# Filter unique samples
setkey(database.info, participant, cancer.barcode, normal.barcode)
database.info.unq <- unique(database.info)

# Assign systematic names to data frames and save ------------------------------
write.table(database.info.unq, paste(parsed.data.dir, cancer.type,
                                     "/info/muta-participants.txt", sep=""),
            quote=F, row.names=F)
