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
mutation <- fread(paste(source.file.dir, file.sample.map[1,filename], sep=""), 
                  select=c(1,2,5,6,7,9,16,17))
setnames(mutation, c("gene", "entrez.id", "chrom", "start", "end", 
                     "type", "cancer.barcode", "normal.barcode"))

# Assign systematic names to data frames and save ------------------------------
assign(paste(cancer.type, ".muta.cancer", sep=""), mutation)
save(list = paste(cancer.type, ".muta.cancer", sep=""),
     file = paste(output.dir, cancer.type, "-muta-cancer.Rdata", sep=""))
