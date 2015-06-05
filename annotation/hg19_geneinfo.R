library(data.table)

cancer_type <- "brca"
raw_data_dir    <- "/share/scratch/arj32/raw_data/"
parsed_data_dir <- "/share/scratch/arj32/parsed_data/"

raw <- fread(paste(raw_data_dir, "annotation/RefSeqHg19.txt", sep=""),
             select=c(3,5,6,13))

filter <- raw[,.(start=min(txStart), end=max(txEnd)),by=.(chrom,name2)]
filter[,chrom:=gsub("^chr","",chrom)]

geneInfo <- filter[chrom %in% c(1:23,"X","Y")]

setnames(geneInfo, c("chrom", "genename", "start", "end"))
setcolorder(geneInfo, c("chrom", "start", "end", "genename"))
setkey(geneInfo, chrom, start, end)

write.table(geneInfo,
            paste(parsed_data_dir, cancer_type,
                  "/annotation/hg19geneinfo.txt", sep=""),
            quote=F, row.names=F)
