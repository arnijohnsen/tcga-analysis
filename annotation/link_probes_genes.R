library(data.table)

cancer_type <- "brca"
raw_data_dir <- "/share/scratch/arj32/raw_data/"
parsed_data_dir <- "/share/scratch/arj32/parsed_data/"

cat("Reading probes annotation file..\n")
probe_ann <- fread(paste(raw_data_dir,
                         "annotation/GenomeStudioProbeAnnotations.txt",
                         sep=""),
                   select=c(2, 12))
setnames(probe_ann, c("probe", "gene"))

gene_list <- strsplit(probe_ann$gene, ";")
probe_list <- rep(probe_ann$probe, sapply(gene_list, length))

probe_ann_split <- unique(data.table(probe = probe_list,
                                     gene  = unlist(gene_list)))

write.table(probe_ann_split, paste(parsed_data_dir, cancer_type,
                                   "/annotation/meth_linked_probes_genes.txt",
                                   sep = ""),
            quote = F, row.names = F)
