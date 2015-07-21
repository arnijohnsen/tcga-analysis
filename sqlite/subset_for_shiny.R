library(data.table)

source_dir <- "/share/scratch/arj32/parsed_data/"
output_dir <- "/share/scratch/arj32/shiny_data/"
cancer_type <- "brca"
data_type <- "lpg"

nice_genes <- unlist(
  fread(
    "/share/scratch/arj32/tcga-analysis/sqlite/nice_genes.txt",
    header = F
  ),
  use.names = F
)
nice_mirns <- unlist(
  fread(
    "/share/scratch/arj32/tcga-analysis/sqlite/nice_mirns.txt",
    header = F
  ),
  use.names = F
)

# Load cnvs and save information for nice_genes one by one ---------------------
if (data_type == "cnvs") {
  cat("Subsetting copy number\n")
  cnvs <- readRDS(
    paste(
      source_dir, cancer_type, "/cnvs/", cancer_type, "_cnvs_cancer.Rds", sep = ""
    )
  )
  setkey(cnvs, gene)

  cnvs_output_dir <- paste(output_dir, cancer_type, "/cnvs/", sep = "")
  for (gene in c(nice_genes, nice_mirns)) {
    tmp <- unlist(cnvs[gene, -(1:4), with = F], use.names = F)
    if (!all(is.na(tmp))) {
      saveRDS(
        tmp, 
        file = paste(cnvs_output_dir, gene, ".Rds", sep = "")
      )
    }
  }
  cnvs_participants <- substr(names(cnvs)[-(1:4)], 1, 12)
  saveRDS(
    cnvs_participants,
    file = paste(output_dir, cancer_type, "/cnvs_participants.Rds", sep = "")
  )
}

# Load expr and save information for nice_genes one by one ---------------------
if (data_type == "expr") {
  cat("Subsetting expression\n")
  expr <- readRDS(
    paste(
      source_dir, cancer_type, "/expr/", cancer_type, "_expr_cancer.Rds", sep = ""
    )
  )
  setkey(expr, gene)

  expr_output_dir <- paste(output_dir, cancer_type, "/expr/", sep = "")
  for (gene in nice_genes) {
    tmp <- unlist(expr[gene, -1, with = F], use.names = F)
    if (!all(is.na(tmp))) {
      saveRDS(
        tmp, 
        file = paste(expr_output_dir, gene, ".Rds", sep = "")
      )
    }
  }
  expr_participants <- substr(names(expr)[-1], 1, 12)
  saveRDS(
    expr_participants,
    file = paste(output_dir, cancer_type, "/expr_participants.Rds", sep = "")
  )
}

# Load mirn and save information for nice_mirns one by one ---------------------
if (data_type == "mirn") {
  cat("Subsetting miRNA expression\n")
  mirn <- readRDS(
    paste(
      source_dir, cancer_type, "/mirn/", cancer_type, "_mirn_cancer.Rds", sep = ""
    )
  )
  mirn[,mirn:=toupper(gsub("hsa-", "", mirn))]
  mirn[,mirn:=gsub("-", "", mirn)]
  setkey(mirn, mirn)

  mirn_output_dir <- paste(output_dir, cancer_type, "/expr/", sep = "")
  for (gene in nice_mirns) {
    tmp <- unlist(mirn[gene, -1, with = F], use.names = F)
    if (!all(is.na(tmp))) {
      saveRDS(
        tmp, 
        file = paste(mirn_output_dir, gene, ".Rds", sep = "")
      )
    }
  }
  mirn_participants <- substr(names(mirn)[-1], 1, 12)
  saveRDS(
    mirn_participants,
    file = paste(output_dir, cancer_type, "/mirn_participants.Rds", sep = "")
  )
}

# Load meth and save information for probes associated to nice_genes -----------
if (data_type == "meth") {
  cat("Subsetting methylation\n")
  meth <- readRDS(
    paste(
      source_dir, cancer_type, "/meth/", cancer_type, "_meth_cancer.Rds",
      sep = ""
    )
  )
  setkey(meth, probe)

  linked_probes_genes <- fread(
    paste(
      source_dir, cancer_type, "/info/meth_linked_probes_genes.txt", sep = ""
    )
  )

  linked_probes_genes[,gene:=gsub("MIRLET", "LET", gene)]
  setkey(linked_probes_genes, gene)

  meth_output_dir <- paste(output_dir, cancer_type, "/meth/", sep = "")
  for (gene in c(nice_genes, nice_mirns)) {
    probes <- linked_probes_genes[gene]$probe
    for (probe in probes) {
      tmp <- unlist(meth[probe, -(1:4), with = F], use.names = F)
      if (!all(is.na(tmp))) {
        saveRDS(
          tmp, 
          file = paste(meth_output_dir, probe, ".Rds", sep = "")
        )
      }
    }
  }
  meth_participants <- substr(names(meth)[-(1:4)], 1, 12)
  saveRDS(
    meth_participants,
    file = paste(output_dir, cancer_type, "/meth_participants.Rds", sep = "")
  )
}

# Load muta and save information for nice_genes one by one ---------------------
if (data_type == "muta") {
  cat("Subsetting mutation\n")
  muta <- readRDS(
    paste(
      source_dir, cancer_type, "/muta/", cancer_type, "_muta_cancer.Rds", sep = ""
    )
  )
  muta[,gene:=gsub("-", "", gene)]
  setkey(muta, gene)

  muta_output_dir <- paste(output_dir, cancer_type, "/muta/", sep = "")
  for (gene in c(nice_genes, nice_mirns)) {
    tmp <- muta[gene]
    tmp2 <- data.frame(
      participant = substr(tmp$cancer_barcode, 1, 12),
      type = tmp$type,
      stringsAsFactors = FALSE
    )
    if (all(is.na(tmp2))) {
      tmp2 <- NULL
    }
    saveRDS(
      tmp2,
      file = paste(muta_output_dir, gene, ".Rds", sep = "")
    )
  }
  muta_participants <- unique(substr(muta$cancer_barcode, 1, 12))
  saveRDS(
    muta_participants,
    file = paste(output_dir, cancer_type, "/muta_participants.Rds", sep = "")
  )
}

if (data_type == "clin") {
  clin <- readRDS(
    paste(
      source_dir, cancer_type, "/clin/", cancer_type, "_clin_cancer.Rds", sep = ""
    )
  )
  tmp <- data.frame(clin, stringsAsFactors = FALSE)
  saveRDS(
    tmp,
    file = paste(output_dir, cancer_type, "/clin.Rds", sep = "")
  )
}

if (data_type == "lpg") {
  linked_probes_genes <- fread(
    paste(
      source_dir, cancer_type, "/info/meth_linked_probes_genes.txt", sep = ""
    )
  )
  linked_probes_genes[,gene:=gsub("MIRLET", "LET", gene)]
  setkey(linked_probes_genes, gene)
  probe_status <- fread(
    paste(
      source_dir, cancer_type, "/info/meth_probe_status.txt", sep = ""
    )
  )
  setkey(probe_status, probe)
  tmp <- data.frame(linked_probes_genes[c(nice_genes, nice_mirns)])
  tmp$status <- probe_status[tmp$probe]$status
  tmp$status[is.na(tmp$status)] <- "undefined"
  saveRDS(
    tmp, 
    file = paste(output_dir, cancer_type, "/linked_probes_genes.Rds", sep="")
  )
}
