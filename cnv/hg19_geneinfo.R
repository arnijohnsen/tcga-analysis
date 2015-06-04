library(data.table)

cancer.type <- "brca"
raw.data.dir    <- "/share/scratch/arj32/raw_data/"
parsed.data.dir <- "/share/scratch/arj32/parsed_data/"

raw <- fread(paste(raw.data.dir, "annotation_files/RefSeqHg19.txt", sep=""),
             select=c(3,5,6,13))

filter <- raw[,.(start=min(txStart), end=max(txEnd)),by=.(chrom,name2)]
filter[,chrom:=gsub("^chr","",chrom)]

geneInfo <- filter[chrom %in% c(1:23,"X","Y")]

setnames(geneInfo, c("chrom", "genename", "start", "end"))
setcolorder(geneInfo, c("chrom", "start", "end", "genename"))
setkey(geneInfo, chrom, start, end)

write.table(geneInfo,
            paste(parsed.data.dir, cancer.type,
                  "/info/hg19geneinfo.txt", sep=""),
            quote=F, row.names=F)
