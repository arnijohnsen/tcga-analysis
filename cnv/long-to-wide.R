library(CNTools)
library(data.table)

data(geneInfo) # based on build 36

load("../parsed-data/OV/cnv/ov-cnv-cancer.Rdata")
cnv.long <- rbindlist(ov.cnv.cancer[1:10])
cnv.long.filter <- cnv.long[Segment_Mean < -0.1 | Segment_Mean > 0.1]
setnames(cnv.long.filter, c("ID", "chrom", "loc.start", "loc.end", "n.probes", "seg.mean"))
#cnv.long.filter[,n.probes:=NULL]
cnv.long.filter$ID <- substr(cnv.long.filter$ID, 6, 12)

cn.seg <- CNSeg(cnv.long.filter)
#rd.seg <- getRS(cn.seg, by="region", imput=F, XY=T, what="mean")
rd.seg <- getRS(cn.seg, by="gene", imput=F, XY=T, what="mean", geneMap=geneInfo)

cnv.wide <-data.table(rs(rd.seg))
#cnv.wide.filtered <- cnv.wide[,.SD[!all(.SD == 0)], by=.(chrom,start,end)]
