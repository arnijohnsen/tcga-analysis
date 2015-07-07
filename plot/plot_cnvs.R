library(copynumber)
library(data.table)

cat("..loading data\n")
brca_cnvl_cancer <- readRDS("../parsed_data/brca/cnvs/brca_cnvl_cancer.Rds")

cat("..binding list\n")
cnvl_bound <- rbindlist(brca_cnvl_cancer[1:20])

cat("..setting names and keys\n")
setnames(
  cnvl_bound,
  c("sampleID", "chrom", "start.pos", "end.pos", "n.probes", "mean")
)
setkey(cnvl_bound, chrom, start.pos, end.pos)

cat("..creating data.frame\n")
x <- data.frame(cnvl_bound)

cat("..plotting\n")
plotFreq(x, thres.gain=0.2, thres.loss=-0.2)
