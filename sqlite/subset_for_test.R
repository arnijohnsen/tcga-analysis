library(data.table)
library(RSQLite)
library(yaml)    # Read config
library(utils)   # Progress bar
library(R.utils) # Verbose

# Define verbose and verbose level ---------------------------------------------
verbose <- Verbose(threshold = -1)

# Read base directories from config.yml, define cancer and data type -----------
ruler(verbose)
enter(verbose, "Preparation")
cat(verbose, "Reading config.yml")
config <- yaml.load_file("/share/scratch/arj32/tcga-analysis/config.yml")
raw_data_dir    <- config$data_dirs$raw_data
parsed_data_dir <- config$data_dirs$parsed_data
cancer_type     <- config$cancer_type

# Define list of interesting genes ---------------------------------------------
cat(verbose, "Defining interesting genes")
nice_genes <- c(
  "UNG",    "SMUG1",  "MBD4",    "TDG",     "OGG1",   "MUTYH",   "NTHL1",
  "MPG",    "NEIL1",  "NEIL2",   "NEIL3",   "APEX1",  "APEX2",   "LIG3",
  "XRCC1",  "PNKP",   "APLF",    "PARP1",   "PARP2",  "PARP3",   "MGMT",
  "ALKBH2", "ALKBH3", "TDP1",    "MSH2",    "MSH3",   "MSH6",    "MLH1",
  "PMS2",   "MSH4",   "MSH5",    "MLH3",    "PMS1",   "PMS2L3",  "XPC",
  "RAD23B", "CETN2",  "RAD23A",  "XPA",     "DDB1",   "DDB2",    "RPA1",
  "RPA2",   "RPA3",   "ERCC3",   "ERCC2",   "GTF2H1", "GTF2H2C", "GTF2H2D",
  "GTF2H3", "GTF2H4", "GTF2H5",  "CDK7",    "CCNH",   "MNAT1",   "ERCC5",
  "ERCC1",  "ERCC4",  "LIG1",    "ERCC8",   "ERCC6",  "XAB2",    "MMS19",
  "RAD51",  "DMC1",   "XRCC2",   "XRCC3",   "RAD52",  "RAD54L",  "RAD54B",
  "BRCA1",  "SHFM1",  "RAD50",   "MRE11A",  "NBN",    "RBBP8",   "MUS81",
  "EME1",   "EME2",   "GIYD1",   "GIYD2",   "GEN1",   "FANCA",   "FANCB",
  "FANCC",  "BRCA2",  "FANCD2",  "FANCE",   "FANCF",  "FANCG",   "FANCI",
  "BRIP1",  "FANCL",  "FANCM",   "PALB2",   "RAD51C", "BTBD12",
  "XRCC6",  "XRCC5",  "PRKDC",   "LIG4",    "XRCC4",  "DCLRE1C", "NHEJ1",
  "NUDT1",  "DUT",    "RRM2B",   "POLB",    "POLG",   "POLD1",   "POLE",
  "PCNA",   "REV3L",  "MAD2L2",  "POLH",    "POLI",   "POLQ",
  "POLK",   "POLL",   "POLM",    "POLN",    "FEN1",   "TREX1",
  "TREX2",  "EXO1",   "APTX",    "SPO11",   "UBE2A",  "UBE2B",
  "RAD18",  "SHPRH",  "HLTF",    "RNF168",  "RNF8",   "RNF4",
  "UBE2V2", "UBE2N",  "H2AFX",   "CHAF1A",  "SETMAR", "BLM",     "WRN",
  "RECQL4", "ATM",    "DCLRE1A", "DCLRE1B", "RPA4",
  "PRPF19", "RECQL",  "RECQL5",  "HELQ",    "RDM1",   "OBFC2B",  "ATR",
  "ATRIP",  "MDC1",   "RAD1",    "RAD9A",   "HUS1",   "RAD17",   "CHEK1",
  "CHEK2",  "TP53",   "TP53BP1", "RIF1",    "TOPBP1", "CLK2",    "PER1",
  "ERBB2",  "CDH1",   "RINT1"
)
exit(verbose)

# Open large data files --------------------------------------------------------
enter(verbose, "Loading .Rds files")
cat(verbose, "Clinical data")
clin <- readRDS(
  paste(parsed_data_dir, cancer_type, "/clin/brca_clin_cancer.Rds", sep = "")
)
cat(verbose, "Copy number data")
cnvs_cancer <- readRDS(
  paste(parsed_data_dir, cancer_type, "/cnvs/brca_cnvs_cancer.Rds", sep = "")
)
cat(verbose, "Expression data")
expr_cancer <- readRDS(
  paste(parsed_data_dir, cancer_type, "/expr/brca_expr_cancer.Rds", sep = "")
)
expr_normal <- readRDS(
  paste(parsed_data_dir, cancer_type, "/expr/brca_expr_normal.Rds", sep = "")
)
cat(verbose, "Methylation data")
meth_cancer <- readRDS(
  paste(parsed_data_dir, cancer_type, "/meth/brca_meth_cancer.Rds", sep = "")
)
meth_normal <- readRDS(
  paste(parsed_data_dir, cancer_type, "/meth/brca_meth_normal.Rds", sep = "")
)
cat(verbose, "Mutation data")
muta_cancer <- readRDS(
  paste(parsed_data_dir, cancer_type, "/muta/brca_muta_cancer.Rds", sep = "")
)
exit(verbose)

# Subset data with nice_genes --------------------------------------------------
enter(verbose, "Subsetting data with nice genes")
cat(verbose, "Clinical is not subsetted")
cat(verbose, "Copy number")
setkey(cnvs_cancer, gene)
cnvs_cancer_min <- cnvs_cancer[nice_genes]
cat(verbose, "Expression")
setkey(expr_cancer, gene)
setkey(expr_normal, gene)
expr_cancer_min <- expr_cancer[nice_genes]
expr_normal_min <- expr_normal[nice_genes]
cat(verbose, "Methylation")
lpg <- fread(
  paste(
    parsed_data_dir, cancer_type, "/info/meth_linked_probes_genes.txt", sep = ""
  )
)
setkey(lpg, gene)
setkey(meth_cancer, probe)
setkey(meth_normal, probe)
nice_probes <- intersect(lpg[nice_genes]$probe, meth_cancer$probe)
meth_cancer_min <- meth_cancer[nice_probes]
meth_normal_min <- meth_normal[nice_probes]
cat(verbose, "Mutation is not subsetted")
exit(verbose)

# Save to database -------------------------------------------------------------
enter(verbose, "Saving to database")
db_name <- paste(
  parsed_data_dir, cancer_type, "/", cancer_type, "_min.sqlite", sep=""
)
cat(verbose, "Opning databse at ", db_name)
db <- dbConnect(SQLite(), dbname = db_name)

if(length(dbListTables(db)) != 0){
  stop("Database isn't empty\n")
}

enter(verbose, "Sending writing queries for")
cat(verbose, "Gene names")
dbWriteTable(
  conn = db, name = "nice_genes", value = data.frame(nice_genes), row.names=F
)
cat(verbose, "Probe names")
dbWriteTable(
  conn = db, name = "nice_probes", value = data.frame(nice_probes), row.names=F
)
dbWriteTable(
  conn = db, name = "linked_probes_genes", value = lpg[nice_genes], row.names=F
)
cat(verbose, "Clinical")
dbWriteTable(
  conn = db, name = "clin", value = clin, row.names = F
)
cat(verbose, "Copy number")
dbWriteTable(
  conn = db, name = "cnvs_cancer", value = cnvs_cancer_min, row.names = F
)
cat(verbose, "Expression")
dbWriteTable(
  conn = db, name = "expr_cancer", value = expr_cancer_min, row.names = F
)
dbWriteTable(
  conn = db, name = "expr_normal", value = expr_normal_min, row.names = F
)
cat(verbose, "Methylation")
dbWriteTable(
  conn = db, name = "meth_cancer", value = meth_cancer_min, row.names = F
)
dbWriteTable(
  conn = db, name = "meth_normal", value = meth_normal_min, row.names = F
)
cat(verbose, "Mutation")
dbWriteTable(
  conn = db, name = "muta_cancer", value = muta_cancer, row.names = F
)
dbDisconnect(db)
