library(data.table)

cancer.type <- "brca"
raw.data.dir <- "/share/scratch/arj32/raw-data/"
parsed.data.dir <- "/share/scratch/arj32/parsed-data/"

cat("Reading probes annotation file..\n")
probe.ann <- fread(paste(raw.data.dir, 
                         "annotation-files/GenomeStudioProbeAnnotations.txt",
                         sep=""),
                    select=c(2,14, 16,19))
setnames(probe.ann, c("probe", "group", "cpg.island", "is.enhancer"))

cat("Setting probe status..\n")
# Creating logical vectors which indicate which probe is in which category -----
probe.ann[,is.body  :=!grepl("TSS200|TSS1500|1stExon|5'UTR|3'UTR|^$", group)]
probe.ann[,is.island:= grepl("Island", cpg.island)]
probe.ann[,is.shore := grepl("Shore",  cpg.island)]
probe.ann[,is.none  :=(!is.island & !is.shore)]
probe.ann[,is.promoter:=grepl("TSS200|5'UTR", group)]
probe.ann[is.na(is.enhancer)]$is.enhancer <- FALSE

# Creating status column, with default case (undefined) and then adding other --
# status types in increasing order of importance (most important is laste) -----
probe.ann[,status:="undefined"]
probe.ann[is.body & is.island]$status <- "body.island"
probe.ann[is.body & is.shore ]$status <- "body.shore"
probe.ann[is.body & is.none  ]$status <- "body.none"
probe.ann[(is.enhancer)      ]$status <- "enhancer"
probe.ann[(is.promoter)      ]$status <- "promoter"

cat("Writing to file..\n")
write.table(probe.ann[,c(1,10),with=F],
            paste(parsed.data.dir,
                  cancer.type,
                  "/info/meth-probe-status.txt",
                  sep=""),
            quote=F, row.names=F)
