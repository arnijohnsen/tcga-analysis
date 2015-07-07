library(data.table)
library(CNTools)

n <- 10

if (!exists("cnvl")) {
  cnvl <- readRDS("../parsed_data/brca/cnvs/brca_cnvl_cancer.Rds")
}

chrom_lengths <- c(
  249250621,
  243199373,
  198022430,
  191154276,
  180915260,
  171115067,
  159138663,
  146364022,
  141213431,
  135534747,
  135006516,
  133851895,
  115169878,
  107349540,
  102531392,
  90354753,
  81195210,
  78077248,
  59128983,
  63025520,
  48129895,
  51304566,
  155270560,
  59373566
)

chrom_pos_to_abs_pos <- function(chrom, pos) {
  chrom[chrom == "X"] <- 23
  chrom[chrom == "Y"] <- 24
  chrom <- as.integer(chrom)
  return(
    sapply(
      chrom,
      function(x) {
        sum(as.numeric(chrom_lengths[1:x-1]))
      }
    ) + pos
  )
}

cnvl_part <- rbindlist(cnvl[1:n])
setnames(cnvl_part, c("ID", "chrom", "loc.start", "loc.end", "n.probes", "seg.mean"))
cnvs_part <- data.table(
  rs(
    getRS(
      CNSeg(cnvl_part),
      by = "region", imput = F, XY = T, what = "mean"
    )
  )
)
cnvs_part[, chrom:=as.character(chrom)]
cnvs_part[, start:=as.numeric(levels(start))[start]]
cnvs_part[, end  :=as.numeric(levels(end)  )[end]]


cnvs_part$start_abs <- chrom_pos_to_abs_pos(cnvs_part$chrom, cnvs_part$start)

plot_colors <- rainbow(n)
x <- unlist(cnvs_part[,4,with=F])
plot(
  cnvs_part$start_abs, as.numeric(levels(x))[x], type="s", col=plot_colors[1],
  xlim = c(7.6e8,7.8e8)
)
abline(v = sapply(1:length(chrom_lengths),function(x) {sum(chrom_lengths[1:x])}), col="gray")
for (i in 5:(n+3)) {
  x <- unlist(cnvs_part[,i,with=F])
  lines(cnvs_part$start_abs, as.numeric(levels(x))[x], type="s",col=plot_colors[i-3])
}
