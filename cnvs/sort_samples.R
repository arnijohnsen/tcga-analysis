library(data.table)

cancer_type     <- "brca"
data_type       <- "cnv"
raw_data_dir    <- "/share/scratch/arj32/raw_data/"
parsed_data_dir <- "/share/scratch/arj32/parsed_data/"

# Read FILE_SAMPLE_MAP.txt and filter correct type of file ---------------------
file_sample_map <- fread(paste(raw_data_dir, cancer_type,
                               "/cnv/FILE_SAMPLE_MAP.txt", sep="")
                        )[grep("nocnv_hg19", filename)]

setnames(file_sample_map, c("filename", "barcode"))

# Remove all samples which aren't primary solid tumor ("01") -------------------
file_sample_map <- file_sample_map[substr(barcode, 14, 15) %in% c("01")]
# file_sample_map[, filename := NULL]
file_sample_map[, participant := substr(barcode, 1, 12)]
setkey(file_sample_map, participant)

file_sample_map <- file_sample_map[!grepl(",", barcode)]
file_sample_map[,filename:=NULL]

setcolorder(file_sample_map, c("participant", "barcode"))
setnames(file_sample_map, c("participant", "cancer_barcode"))


# multi_barcode <- file_sample_map[grep(",", barcode)]
# multi_barcode[, name:= make.names(participant, unique=T)]

# x <- list()

# for(i in 1:dim(multi_barcode)[1]){
#   x[[multi_barcode$name[i]]] <- fread(paste(raw_data_dir, "brca/cnv/CNV_SNP_Array/BI__Genome_Wide_SNP_6/Level_3/", multi_barcode$filename[i], sep=""))
#   x[[i]]$Sample <- multi_barcode$name[i]
# }
# w <- rbindlist(x)

# setnames(w, c("sampleID", "chrom", "start.pos", "end.pos", "n.probes", "mean"))
# w[chrom=="X"]$chrom <- "23"
# w[,chrom:=as.numeric(chrom)]

# setkey(w, sampleID, chrom, start.pos, end.pos)
# w[, arm:=rep("p", dim(w)[1])]
# setcolorder(w, c("sampleID", "chrom", "arm", "start.pos", "end.pos", "n.probes", "mean"))

# w_df <- data.frame(w)
# plotHeatmap(segments=w_df, upper.lim=0.1)
