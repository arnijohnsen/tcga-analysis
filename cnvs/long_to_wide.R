library(CNTools)
library(data.table)

cancer_type <- "brca"
raw_data_dir    <- "/share/scratch/arj32/raw_data/"
parsed_data_dir <- "/share/scratch/arj32/parsed_data/"
output_dir      <- paste(parsed_data_dir, cancer_type, "/cnvs/",  sep="")

#data(geneInfo) # based on build 36
# Load hg19 (build 37)
geneInfo <- fread(paste(parsed_data_dir, cancer_type, "/info/hg19geneinfo.txt", sep=""))

load(paste(output_dir, cancer_type, "_cnvl_cancer.Rdata", sep=""))
load(paste(output_dir, cancer_type, "_cnvl_normal.Rdata", sep=""))

cat("Binding list..\n")
cnvl_bound_cancer <- rbindlist(brca_cnvl_cancer)
cnvs_bound_normal <- rbindlist(brca_cnvl_normal)

setnames(cnvl_bound_cancer, c("ID", "chrom", "loc.start", "loc.end", "n.probes", "seg.mean"))
setnames(cnvs_bound_normal, c("ID", "chrom", "loc.start", "loc.end", "n.probes", "seg.mean"))

cat("CNSeg..\n")
cn_seg_cancer <- CNSeg(cnvl_bound_cancer)
cn_seg_normal <- CNSeg(cnvs_bound_normal)

cat("RDSeg..\n")
rd_seg_cancer <- getRS(cn_seg_cancer, by="gene", imput=F, XY=T, what="mean", geneMap=geneInfo)
rd_seg_normal <- getRS(cn_seg_normal, by="gene", imput=F, XY=T, what="mean", geneMap=geneInfo)

cat("rs(RDSeg)..\n")
cnvs_cancer <-data.table(rs(rd_seg_cancer))
cnvs_normal <-data.table(rs(rd_seg_normal))

assign(paste(cancer_type, "_cnvs_cancer", sep=""), cnvs_cancer)
assign(paste(cancer_type, "_cnvs_normal", sep=""), cnvs_normal)
save(list = paste(cancer_type, "_cnvs_cancer", sep=""),
     file = paste(output_dir, cancer_type, "_cnvs_cancer.Rdata", sep=""))
save(list = paste(cancer_type, "_cnvs_normal", sep=""),
     file = paste(output_dir, cancer_type, "_cnvs_normal.Rdata", sep=""))
