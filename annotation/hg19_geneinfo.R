library(data.table)

cancer_type <- "brca"
raw_data_dir    <- "/share/scratch/arj32/raw_data/"
parsed_data_dir <- "/share/scratch/arj32/parsed_data/"

raw_gene <- fread(
  paste(raw_data_dir, "annotation/RefSeqHg19.txt", sep=""),
  select=c(3,5,6,13)
)
raw_mirn <- fread(
  paste(raw_data_dir, "annotation/hsa.gff3.txt", sep=""),
  select=c(1,3,4,5,9)
)
setnames(raw_mirn, c("chrom", "type", "txStart", "txEnd", "genename"))

raw_mirn <- raw_mirn[type=="miRNA_primary_transcript"]
raw_mirn[,type:=NULL]
raw_mirn[,genename:=gsub('.*Name=hsa-', "", genename)]
raw_mirn[,genename:=gsub('";', "", genename)]
raw_mirn[,genename:=gsub("-", "", genename)]
raw_mirn[,genename:=toupper( genename)]
raw_mirn[,chrom:=gsub("^chr", "", chrom)]

# Remove duplicate of MIR511
filter_mirn <- raw_mirn[,.(start=min(txStart), end=max(txEnd)),by=.(chrom,genename)]
filter_gene <- raw_gene[,.(start=min(txStart), end=max(txEnd)),by=.(chrom,name2)]
filter_gene[,chrom:=gsub("^chr","",chrom)]

gene_info <- filter_gene[chrom %in% c(1:23,"X","Y")]

setnames(gene_info, c("chrom", "genename", "start", "end"))
setcolorder(gene_info, c("chrom", "start", "end", "genename"))

gene_info <- rbind(gene_info, filter_mirn)
setkey(gene_info, chrom, start, end)

write.table(
  gene_info,
  paste(
    parsed_data_dir, cancer_type,
    "/info/hg19geneinfo.txt", sep=""
  ),
  quote=F, row.names=F
)
