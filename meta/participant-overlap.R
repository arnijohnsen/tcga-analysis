library(data.table)

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


cnv.c  <-  cnv[!is.na(cancer.barcode)]
expr.c <- expr[!is.na(cancer.barcode)]
meth.c <- meth[!is.na(cancer.barcode)]

cnv.n  <-  cnv[!is.na(normal.barcode)]
expr.n <- expr[!is.na(normal.barcode)]
meth.n <- meth[!is.na(normal.barcode)]

cnv.b  <-  cnv[!is.na(cancer.barcode)&!is.na(normal.barcode)]
expr.b <- expr[!is.na(cancer.barcode)&!is.na(normal.barcode)]
meth.b <- meth[!is.na(cancer.barcode)&!is.na(normal.barcode)]

overlap <- data.table(c("cnv (total)", "expr (total)", "meth (total)", 
                        "cnv+expr", "cnv+meth", "meth+expr", "all"))
setnames(overlap, c("sets"))

overlap[,total:=c( length( cnv[,participant]),
                   length(expr[,participant]),
                   length(meth[,participant]),
                   length( intersect(cnv[,participant],expr[,participant])),
                   length( intersect(cnv[,participant],meth[,participant])),
                   length(intersect(expr[,participant],meth[,participant])),
                   length(intersect(intersect(cnv[,participant],
                                              expr[,participant]),
                                    meth[,participant])))]
overlap[,cancer:=c( length( cnv.c[,participant]),
                    length(expr.c[,participant]),
                    length(meth.c[,participant]),
                    length( intersect(cnv.c[,participant],expr.c[,participant])),
                    length( intersect(cnv.c[,participant],meth.c[,participant])),
                    length(intersect(expr.c[,participant],meth.c[,participant])),
                    length(intersect(intersect(cnv.c[,participant],
                                               expr.c[,participant]),
                                     meth.c[,participant])))]
overlap[,normal:=c( length( cnv.n[,participant]),
                    length(expr.n[,participant]),
                    length(meth.n[,participant]),
                    length( intersect(cnv.n[,participant],expr.n[,participant])),
                    length( intersect(cnv.n[,participant],meth.n[,participant])),
                    length(intersect(expr.n[,participant],meth.n[,participant])),
                    length(intersect(intersect(cnv.n[,participant],
                                               expr.n[,participant]),
                                     meth.n[,participant])))]
overlap[,  both:=c( length( cnv.b[,participant]),
                    length(expr.b[,participant]),
                    length(meth.b[,participant]),
                    length( intersect(cnv.b[,participant],expr.b[,participant])),
                    length( intersect(cnv.b[,participant],meth.b[,participant])),
                    length(intersect(expr.b[,participant],meth.b[,participant])),
                    length(intersect(intersect(cnv.b[,participant],
                                               expr.b[,participant]),
                                     meth.b[,participant])))]
write.table(overlap, paste(parsed.data.dir, cancer.type, 
            "/info/participant-overlap.txt", sep=""), quote=F, row.names=F)
