library(data.table)
library(limma)

# Define cancer type, raw and parsed data directories --------------------------
cancer.type <- "brca"
raw.data.dir    <- "/share/scratch/arj32/raw-data/"
parsed.data.dir <- "/share/scratch/arj32/parsed-data/"

# Read pariticipant list for cnv, expr and meth. TODO: add muta ----------------
cnv <- fread(paste(parsed.data.dir, cancer.type,
                           "/info/cnv-participants.txt", sep=""))
expr <- fread(paste(parsed.data.dir, cancer.type,
                           "/info/expr-participants.txt", sep=""))
meth <- fread(paste(parsed.data.dir, cancer.type,
                           "/info/meth-participants.txt", sep=""))
muta <- fread(paste(parsed.data.dir, cancer.type,
                           "/info/muta-participants.txt", sep=""))
mirn <- fread(paste(parsed.data.dir, cancer.type,
                           "/info/mirn-participants.txt", sep=""))

cnv.c  <-  cnv[!is.na(cancer.barcode)]
expr.c <- expr[!is.na(cancer.barcode)]
meth.c <- meth[!is.na(cancer.barcode)]
muta.c <- muta[!is.na(cancer.barcode)]
mirn.c <- mirn[!is.na(cancer.barcode)]

cnv.n  <-  cnv[!is.na(normal.barcode)]
expr.n <- expr[!is.na(normal.barcode)]
meth.n <- meth[!is.na(normal.barcode)]
muta.n <- muta[!is.na(normal.barcode)]
mirn.n <- mirn[!is.na(normal.barcode)]

cnv.b  <-  cnv[!is.na(cancer.barcode)&!is.na(normal.barcode)]
expr.b <- expr[!is.na(cancer.barcode)&!is.na(normal.barcode)]
meth.b <- meth[!is.na(cancer.barcode)&!is.na(normal.barcode)]
muta.b <- muta[!is.na(cancer.barcode)&!is.na(normal.barcode)]
mirn.b <- mirn[!is.na(cancer.barcode)&!is.na(normal.barcode)]

all.part <- unique(c(cnv$participant,
                     expr$participant,
                     meth$participant,
                     muta$participant,
                     mirn$participant))

DF.c <- data.table(cnv  = all.part %in% cnv.c$participant,
                   expr = all.part %in% expr.c$participant,
                   meth = all.part %in% meth.c$participant,
                   muta = all.part %in% muta.c$participant,
                   mirn = all.part %in% mirn.c$participant)

DF.n <- data.table(cnv  = all.part %in% cnv.n$participant,
                   expr = all.part %in% expr.n$participant,
                   meth = all.part %in% meth.n$participant,
                   muta = all.part %in% muta.n$participant,
                   mirn = all.part %in% mirn.n$participant)

DF.b <- data.table(cnv  = all.part %in% cnv.b$participant,
                   expr = all.part %in% expr.b$participant,
                   meth = all.part %in% meth.b$participant,
                   muta = all.part %in% muta.b$participant,
                   mirn = all.part %in% mirn.b$participant)

VC.c <- vennCounts(DF.c)
VC.n <- vennCounts(DF.n)
VC.b <- vennCounts(DF.b)
overlap <- data.table(cbind(VC.c, VC.n[,6], VC.b[,6]))
setnames(overlap, c("cnv", "expr", "meth", "muta", "mirn", "c", "n", "b"))

n <- dim(overlap)[1]

overlap$c.sum <- rep(0, n)
overlap$n.sum <- rep(0, n)
overlap$b.sum <- rep(0, n)

for (i in 1:n){
  idx <- sapply(1:32, function(x){
    return(all(overlap[x,1:5,with=F] >= overlap[i,1:5,with=F]))
  })
  overlap$c.sum[i] <- sum(overlap$c[idx])
  overlap$n.sum[i] <- sum(overlap$n[idx])
  overlap$b.sum[i] <- sum(overlap$b[idx])
}

#sink(file=paste(parsed.data.dir, cancer.type, 
#                "/info/participant-overlap.txt", sep=""))
#print(overlap)
#sink()
