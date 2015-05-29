library(copynumber)
library(data.table)

cat("..loading data\n")
load("../parsed-data/OV/cnv/ov-cnv-cancer.Rata")

cat("..binding list\n")
cnv.long <- rbindlist(ov.cnv.cancer[1:20])

cat("..setting names and keys\n")
setnames(cnv.long, c("sampleID", "chrom", "start.pos", "end.pos", "n.probes", "mean"))
setkey(cnv.long, chrom, start.pos, end.pos)

cat("..creating data.frame\n")
x <- data.frame(cnv.long)

cat("plotting\n")
plotFreq(x, thres.gain=0.2, thres.loss=-0.2)
