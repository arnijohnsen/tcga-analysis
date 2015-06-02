library(data.table)

# Define cancer type, raw and parsed data directories --------------------------
cancer.type <- "brca"
raw.data.dir    <- "/share/scratch/arj32/raw-data/"
parsed.data.dir <- "/share/scratch/arj32/parsed-data/"
source.file.path<- "/cnv/CNV_SNP_Array/BI__Genome_Wide_SNP_6/Level_3/"
source.file.dir <- paste(raw.data.dir, cancer.type, source.file.path, sep="")
output.dir      <- paste(parsed.data.dir, cancer.type, "/cnv/",  sep="")

# Load file and sample information ---------------------------------------------
file.sample.map <- fread(paste(raw.data.dir, cancer.type, 
                               "/cnv/FILE_SAMPLE_MAP.txt", sep=""))
setnames(file.sample.map, c("filename", "barcode"))
# Use only nocsv_hg19
file.sample.map <- file.sample.map[grep("nocnv_hg19", filename)]

sample.list <- fread(paste(parsed.data.dir, cancer.type, 
                           "/info/cnv-participants.txt", sep=""))

n <- dim(sample.list)[1]

# Create empty lists -----------------------------------------------------------
cancer.list <- list()
normal.list <- list()

# Read files and create a list of data.tables ----------------------------------
for(i in 1:n){
  cat("Reading files for participant", i, "of", n, "\n")
  cancer.barcode <- sample.list$cancer.barcode[i]
  normal.barcode <- sample.list$normal.barcode[i]
  cancer.file <- file.sample.map[barcode==cancer.barcode, filename]
  normal.file <- file.sample.map[barcode==normal.barcode, filename]
  if(!is.na(cancer.barcode)){
    tmp <- fread(paste(source.file.dir, cancer.file, sep=""))
    tmp$Sample <- cancer.barcode
    cancer.list[[ cancer.barcode ]] <- tmp
  }
  if(!is.na(normal.barcode)){
    tmp <- fread(paste(source.file.dir, normal.file, sep=""))
    tmp$Sample <- normal.barcode
    normal.list[[ normal.barcode]] <- tmp
  }
}

# Assign systematic names to data frames and save ------------------------------
assign(paste(cancer.type, ".cnvl.cancer", sep=""), cancer.list)
assign(paste(cancer.type, ".cnvl.normal", sep=""), normal.list)
save(list = paste(cancer.type, ".cnvl.cancer", sep=""),
     file = paste(output.dir, cancer.type, "-cnvl-cancer.Rdata", sep=""))
save(list = paste(cancer.type, ".cnvl.normal", sep=""),
     file = paste(output.dir, cancer.type, "-cnvl-normal.Rdata", sep=""))
