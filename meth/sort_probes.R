library(data.table)

cancer_type <- "brca"
raw_data_dir <- "/share/scratch/arj32/raw_data/"
parsed_data_dir <- "/share/scratch/arj32/parsed_data/"

cat("Reading probes annotation file..\n")
probe_ann <- fread(paste(raw_data_dir,
                         "annotation_files/GenomeStudioProbeAnnotations.txt",
                         sep=""),
                    select=c(2,14, 16, 19))
setnames(probe_ann, c("probe", "group", "cpg_island", "is_enhancer"))

cat("Setting probe status..\n")
# Creating logical vectors which indicate which probe is in which category -----
probe_ann[,is_body  :=!grepl("TSS200|TSS1500|1stExon|5'UTR|3'UTR|^$", group)]
probe_ann[,is_island:= grepl("Island", cpg_island)]
probe_ann[,is_shore := grepl("Shore",  cpg_island)]
probe_ann[,is_none  :=(!is_island & !is_shore)]
probe_ann[,is_promoter:=grepl("TSS200|5'UTR", group)]
probe_ann[is_na(is_enhancer)]$is_enhancer <- FALSE

# Creating status column, with default case (undefined) and then adding other --
# status types in increasing order of importance (most important is laste) -----
probe_ann[,status:="undefined"]
probe_ann[is_body & is_island]$status <- "body_island"
probe_ann[is_body & is_shore ]$status <- "body_shore"
probe_ann[is_body & is_none  ]$status <- "body_none"
probe_ann[(is_enhancer)      ]$status <- "enhancer"
probe_ann[(is_promoter)      ]$status <- "promoter"

cat("Writing to file..\n")
write_table(probe_ann[,c(1,10),with=F],
            paste(parsed_data_dir,
                  cancer_type,
                  "/info/meth_probe_status.txt",
                  sep=""),
            quote=F, row.names=F)
