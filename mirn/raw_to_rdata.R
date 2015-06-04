library(data.table)

# Define cancer type, raw and parsed data directories --------------------------
cancer.type <- "brca"
raw.data.dir    <- "/share/scratch/arj32/raw_data/"
parsed.data.dir <- "/share/scratch/arj32/parsed_data/"
source.file.path<- "/mirn/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3/"
source.file.dir <- paste(raw.data.dir, cancer.type, source.file.path, sep="")
output.dir      <- paste(parsed.data.dir, cancer.type, "/mirn/",  sep="")

# Load file and sample information ---------------------------------------------
file.sample.map <- fread(paste(raw.data.dir, cancer.type,
                               "/mirn/FILE_SAMPLE_MAP.txt", sep=""))
setnames(file.sample.map, c("filename", "barcode"))
# use only genes.normalized_results
file.sample.map <- file.sample.map[grep("mirna.quantification", filename)]

sample.list <- fread(paste(parsed.data.dir, cancer.type,
                               "/info/mirn_participants.txt", sep=""))

sample.list[normal.barcode == ""]$normal.barcode <- NA

n <- dim(sample.list)[1]

# Read one file with all information, use to create data.tables ----------------
cat("Initiating tables..\n")
idx.cancer <- which(!is.na(sample.list$cancer.barcode))[1]
idx.normal <- which(!is.na(sample.list$normal.barcode))[1]
cancer.barcode <- sample.list$cancer.barcode[idx.cancer]
normal.barcode <- sample.list$normal.barcode[idx.normal]
cancer.file <- file.sample.map[barcode==cancer.barcode, filename]
normal.file <- file.sample.map[barcode==normal.barcode, filename]

cat(paste(source.file.dir, cancer.file, sep=""), "\n")
cat(paste(source.file.dir, normal.file, sep=""), "\n")

cancer.wide <- fread(paste(source.file.dir, cancer.file, sep=""),select=1)
normal.wide <- fread(paste(source.file.dir, normal.file, sep=""),select=1)

setnames(cancer.wide, c("mirna"))
setnames(normal.wide, c("mirna"))

# Read all files, drop everything except beta and append as a new column--------
for(i in 1:n){
  cat("Reading files for participant", i, "of", n, "\n")
  cancer.barcode <- sample.list$cancer.barcode[i]
  normal.barcode <- sample.list$normal.barcode[i]
  cancer.file <- file.sample.map[barcode==cancer.barcode, filename]
  normal.file <- file.sample.map[barcode==normal.barcode, filename]
  if(!is.na(cancer.barcode)){
    tmp <- fread(paste(source.file.dir, cancer.file, sep=""), select=3)
    cancer.wide[,c(cancer.barcode):=tmp]
  }
  if(!is.na(normal.barcode)){
    tmp <- fread(paste(source.file.dir, normal.file, sep=""), select=3)
    normal.wide[,c(normal.barcode):=tmp]
  }
}

# Assign systematic names to data frames and save ------------------------------
assign(paste(cancer.type, ".mirn.cancer", sep=""), cancer.wide)
assign(paste(cancer.type, ".mirn.normal", sep=""), normal.wide)
save(list = paste(cancer.type, ".mirn.cancer", sep=""),
     file = paste(output.dir, cancer.type, "_mirn_cancer.Rdata", sep=""))
save(list = paste(cancer.type, ".mirn.normal", sep=""),
     file = paste(output.dir, cancer.type, "_mirn_normal.Rdata", sep=""))
