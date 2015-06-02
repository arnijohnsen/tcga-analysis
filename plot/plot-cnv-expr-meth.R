library(data.table)

if(!exists("brca.cnvw.cancer")){
  load("../parsed-data/brca/cnv/brca-cnvw-cancer.Rdata")
}
if(!exists("brca.expr.cancer")){
  load("../parsed-data/brca/expr/brca-expr-cancer.Rdata")
}
if(!exists("brca.meth.cancer")){
  load("../parsed-data/brca/meth/brca-meth-cancer.Rdata")
}

# TODO: Replace with fread from text file
load("../Rdata/BRCA/info/BRCA-linked-probes-genes.Rdata")

str <- "^ERBB2$"

probes <- (BRCA.linked.probes.genes[grep(str,BRCA.linked.probes.genes$genes),1])

my.probe <- probes[1]

x <- brca.cnvw.cancer[grep(str,genename)]
y <- brca.expr.cancer[grep(str,gene)]
z <- brca.meth.cancer[probe == my.probe]

xv <- unlist(x[,-(1:4),with=F])
yv <- unlist(y[,-1    ,with=F])
zv <- unlist(z[,-(1:4),with=F])

names(xv) <- substr(names(xv), 1, 12)
names(yv) <- substr(names(yv), 1, 12)
names(zv) <- substr(names(zv), 1, 12)
int.sam <- intersect(names(xv),
                     intersect(names(yv),names(zv)))

df <- data.frame(cnv = xv[int.sam], expr = yv[int.sam], meth = zv[int.sam])
plot(df)
