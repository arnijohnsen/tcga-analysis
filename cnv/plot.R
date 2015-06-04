library(copynumber)
library(data.table)

cat("..loading data\n")
load("../parsed_data/OV/cnv/ov_cnv_cancer.Rata")

cat("..binding list\n")
cnv_long <- rbindlist(ov_cnv_cancer[1:20])

cat("..setting names and keys\n")
setnames(cnv_long, c("sampleID", "chrom", "start.pos", "end.pos", "n.probes", "mean"))
setkey(cnv_long, chrom, start.pos, end.pos)

cat("..creating data.frame\n")
x <- data.frame(cnv_long)

cat("plotting\n")
plotFreq(x, thres.gain=0.2, thres.loss=-0.2)
