library(CNTools)
library(data.table)

cancer.type <- "brca"
raw.data.dir    <- "/share/scratch/arj32/raw-data/"
parsed.data.dir <- "/share/scratch/arj32/parsed-data/"
output.dir      <- paste(parsed.data.dir, cancer.type, "/cnv/",  sep="")

#data(geneInfo) # based on build 36
# Load hg19 (build 37)
geneInfo <- fread(paste(parsed.data.dir, cancer.type, "/info/hg19geneinfo.txt", sep=""))

load(paste(output.dir, cancer.type, "-cnv-cancer.Rdata", sep=""))
load(paste(output.dir, cancer.type, "-cnv-normal.Rdata", sep=""))

cat("Binding list..\n")
cnv.long.cancer <- rbindlist(brca.cnv.cancer)
cnv.long.normal <- rbindlist(brca.cnv.normal)

setnames(cnv.long.cancer, c("ID", "chrom", "loc.start", "loc.end", "n.probes", "seg.mean"))
setnames(cnv.long.normal, c("ID", "chrom", "loc.start", "loc.end", "n.probes", "seg.mean"))

cat("CNSeg..\n")
cn.seg.cancer <- CNSeg(cnv.long.cancer)
cn.seg.normal <- CNSeg(cnv.long.normal)

cat("RDSeg..\n")
rd.seg.cancer <- getRS(cn.seg.cancer, by="gene", imput=F, XY=T, what="mean", geneMap=geneInfo)
rd.seg.normal <- getRS(cn.seg.normal, by="gene", imput=F, XY=T, what="mean", geneMap=geneInfo)

cat("rs(RDSeg)..\n")
cnv.wide.cancer <-data.table(rs(rd.seg.cancer))
cnv.wide.normal <-data.table(rs(rd.seg.normal))

assign(paste(cancer.type, ".cnvw.cancer", sep=""), cnv.wide.cancer)
assign(paste(cancer.type, ".cnvw.normal", sep=""), cnv.wide.normal)
save(list = paste(cancer.type, ".cnvw.cancer", sep=""),
     file = paste(output.dir, cancer.type, "-cnvw-cancer.Rdata", sep=""))
save(list = paste(cancer.type, ".cnvw.normal", sep=""),
     file = paste(output.dir, cancer.type, "-cnvw-normal.Rdata", sep=""))
