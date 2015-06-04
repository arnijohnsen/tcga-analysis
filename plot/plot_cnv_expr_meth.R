library(data.table)

if(!exists("brca_cnvw_cancer")){
  load("../parsed_data/brca/cnv/brca_cnvw_cancer.Rdata")
}
if(!exists("brca_expr_cancer")){
  load("../parsed_data/brca/expr/brca_expr_cancer.Rdata")
}
if(!exists("brca_meth_cancer")){
  load("../parsed_data/brca/meth/brca_meth_cancer.Rdata")
}

# TODO: Replace with fread from text file
load("../Rdata/BRCA/info/BRCA-linked-probes-genes.Rdata")

str <- "^ERBB2$"

probes <- (BRCA.linked.probes.genes[grep(str,BRCA.linked.probes.genes$genes),1])

my_probe <- probes[1]

x <- brca_cnvw_cancer[grep(str,genename)]
y <- brca_expr_cancer[grep(str,gene)]
z <- brca_meth_cancer[probe == my_probe]

xv <- unlist(x[,-(1:4),with=F])
yv <- unlist(y[,-1    ,with=F])
zv <- unlist(z[,-(1:4),with=F])

names(xv) <- substr(names(xv), 1, 12)
names(yv) <- substr(names(yv), 1, 12)
names(zv) <- substr(names(zv), 1, 12)
int_sam <- intersect(names(xv),
                     intersect(names(yv),names(zv)))

df <- data.frame(cnv = xv[int_sam], expr = yv[int_sam], meth = zv[int_sam])
plot(df)
