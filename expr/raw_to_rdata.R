library(data.table)

# Define cancer type, raw and parsed data directories --------------------------
cancer.type <- "brca"
raw.data.dir    <- "/share/scratch/arj32/raw-data/"
parsed.data.dir <- "/share/scratch/arj32/parsed-data/"
source.file.path<- "/expr/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/"
source.file.dir <- paste(raw.data.dir, cancer.type, source.file.path, sep="")
output.dir      <- paste(parsed.data.dir, cancer.type, "/expr/",  sep="")

# Load file and sample information ---------------------------------------------
file.sample.map <- fread(paste(raw.data.dir, cancer.type,
                               "/expr/FILE_SAMPLE_MAP.txt", sep=""))
setnames(file.sample.map, c("filename", "barcode"))
# use only genes.normalized_results
file.sample.map <- file.sample.map[grep("genes.normalized_results", filename)]

sample.list <- fread(paste(parsed.data.dir, cancer.type,
                               "/info/expr-participants.txt", sep=""))

n <- dim(sample.list)[1]

# Read one file with all information, use to create data.tables ----------------
cat("Initiating tables..\n")
idx.cancer <- which(!is.na(sample.list$cancer.barcode))[1]
idx.normal <- which(!is.na(sample.list$normal.barcode))[1]
cancer.barcode <- sample.list$cancer.barcode[idx.cancer]
normal.barcode <- sample.list$normal.barcode[idx.normal]
cancer.file <- file.sample.map[barcode==cancer.barcode, filename]
normal.file <- file.sample.map[barcode==normal.barcode, filename]

cancer.wide <- fread(paste(source.file.dir, cancer.file, sep=""))
normal.wide <- fread(paste(source.file.dir, normal.file, sep=""))

cancer.wide[,normalized_count:=NULL]
normal.wide[,normalized_count:=NULL]

setnames(cancer.wide, c("gene"))
setnames(normal.wide, c("gene"))

cancer.wide[,gene:=gsub("\\|.*", "",gene)]
normal.wide[,gene:=gsub("\\|.*", "",gene)]

# Read all files, drop everything except beta and append as a new column--------
for(i in 1:n){
  cat("Reading files for participant", i, "of", n, "\n")
  cancer.barcode <- sample.list$cancer.barcode[i]
  normal.barcode <- sample.list$normal.barcode[i]
  cancer.file <- file.sample.map[barcode==cancer.barcode, filename]
  normal.file <- file.sample.map[barcode==normal.barcode, filename]
  if(!is.na(cancer.barcode)){
    tmp <- fread(paste(source.file.dir, cancer.file, sep=""), select=2)
    cancer.wide[,c(cancer.barcode):=tmp]
  }
  if(!is.na(normal.barcode)){
    tmp <- fread(paste(source.file.dir, normal.file, sep=""), select=2)
    normal.wide[,c(normal.barcode):=tmp]
  }
}

# Assign systematic names to data frames and save ------------------------------
assign(paste(cancer.type, ".expr.cancer", sep=""), cancer.wide)
assign(paste(cancer.type, ".expr.normal", sep=""), normal.wide)
save(list = paste(cancer.type, ".expr.cancer", sep=""),
     file = paste(output.dir, cancer.type, "-expr-cancer.Rdata", sep=""))
save(list = paste(cancer.type, ".expr.normal", sep=""),
     file = paste(output.dir, cancer.type, "-expr-normal.Rdata", sep=""))
