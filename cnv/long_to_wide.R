library(CNTools)
library(data.table)

cancer_type <- "brca"
raw_data_dir    <- "/share/scratch/arj32/raw_data/"
parsed_data_dir <- "/share/scratch/arj32/parsed_data/"
output_dir      <- paste(parsed_data_dir, cancer_type, "/cnv/",  sep="")

#data(geneInfo) # based on build 36
# Load hg19 (build 37)
geneInfo <- fread(paste(parsed_data_dir, cancer_type, "/info/hg19geneinfo.txt", sep=""))

load(paste(output_dir, cancer_type, "_cnv_cancer.Rdata", sep=""))
load(paste(output_dir, cancer_type, "_cnv_normal.Rdata", sep=""))

cat("Binding list..\n")
cnv_long_cancer <- rbindlist(brca_cnv_cancer)
cnv_long_normal <- rbindlist(brca_cnv_normal)

setnames(cnv_long_cancer, c("ID", "chrom", "loc.start", "loc.end", "n.probes", "seg.mean"))
setnames(cnv_long_normal, c("ID", "chrom", "loc.start", "loc.end", "n.probes", "seg.mean"))

cat("CNSeg..\n")
cn_seg_cancer <- CNSeg(cnv_long_cancer)
cn_seg_normal <- CNSeg(cnv_long_normal)

cat("RDSeg..\n")
rd_seg_cancer <- getRS(cn_seg_cancer, by="gene", imput=F, XY=T, what="mean", geneMap=geneInfo)
rd_seg_normal <- getRS(cn_seg_normal, by="gene", imput=F, XY=T, what="mean", geneMap=geneInfo)

cat("rs(RDSeg)..\n")
cnv_wide_cancer <-data_table(rs(rd_seg_cancer))
cnv_wide_normal <-data_table(rs(rd_seg_normal))

assign(paste(cancer_type, ".cnvw_cancer", sep=""), cnv_wide_cancer)
assign(paste(cancer_type, ".cnvw_normal", sep=""), cnv_wide_normal)
save(list = paste(cancer_type, ".cnvw_cancer", sep=""),
     file = paste(output_dir, cancer_type, "_cnvw_cancer.Rdata", sep=""))
save(list = paste(cancer_type, ".cnvw_normal", sep=""),
     file = paste(output_dir, cancer_type, "_cnvw_normal.Rdata", sep=""))
