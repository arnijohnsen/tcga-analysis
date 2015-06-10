library(data.table)
library(RSQLite)

cancer_type <- "brca"
raw_data_dir    <- "/share/scratch/arj32/raw_data/"
parsed_data_dir <- "/share/scratch/arj32/parsed_data/"

# Define list of interesting genes ---------------------------------------------

nice_genes <- c("UNG", "SMUG1", "MBD4", "TDG", "OGG1", "MUTYH", "NTHL1", "MPG",
                "NEIL1", "NEIL2", "NEIL3", "APEX1", "APEX2", "LIG3", "XRCC1",
                "PNKP", "APLF", "PARP1", "PARP2", "PARP3", "MGMT", "ALKBH2",
                "ALKBH3", "TDP1", "MSH2", "MSH3", "MSH6", "MLH1",
                "PMS2", "MSH4", "MSH5", "MLH3", "PMS1", "PMS2L3", "XPC",
                "RAD23B", "CETN2", "RAD23A", "XPA", "DDB1", "DDB2", "RPA1",
                "RPA2", "RPA3", "ERCC3", "ERCC2", "GTF2H1", "GTF2H2C",
                "GTF2H2D", "GTF2H3", "GTF2H4", "GTF2H5", "CDK7", "CCNH",
                "MNAT1", "ERCC5", "ERCC1", "ERCC4", "LIG1", "ERCC8", "ERCC6",
                "XAB2", "MMS19", "RAD51", "DMC1",
                "XRCC2", "XRCC3", "RAD52", "RAD54L", "RAD54B", "BRCA1", "SHFM1",
                "RAD50", "MRE11A", "NBN", "RBBP8", "MUS81", "EME1", "EME2",
                "GIYD1", "GIYD2", "GEN1", "FANCA", "FANCB", "FANCC", "BRCA2",
                "FANCD2", "FANCE", "FANCF", "FANCG", "FANCI", "BRIP1", "FANCL",
                "FANCM", "PALB2", "RAD51C", "BTBD12",
                "XRCC6", "XRCC5", "PRKDC", "LIG4", "XRCC4", "DCLRE1C", "NHEJ1",
                "NUDT1", "DUT", "RRM2B", "POLB", "POLG", "POLD1", "POLE",
                "PCNA", "REV3L", "MAD2L2", "POLH", "POLI", "POLQ",
                "POLK", "POLL", "POLM", "POLN", "FEN1", "TREX1",
                "TREX2", "EXO1", "APTX", "SPO11", "UBE2A", "UBE2B",
                "RAD18", "SHPRH", "HLTF", "RNF168", "RNF8", "RNF4",
                "UBE2V2", "UBE2N", "H2AFX", "CHAF1A", "SETMAR", "BLM", "WRN",
                "RECQL4", "ATM", "DCLRE1A", "DCLRE1B", "RPA4",
                "PRPF19", "RECQL", "RECQL5", "HELQ", "RDM1", "OBFC2B", "ATR",
                "ATRIP", "MDC1", "RAD1", "RAD9A", "HUS1", "RAD17", "CHEK1",
                "CHEK2", "TP53", "TP53BP1", "RIF1", "TOPBP1", "CLK2", "PER1",
                "ERBB2", "CDH1", "RINT1")

# Open large data files --------------------------------------------------------
cat("Loading data..\n")
load(paste(parsed_data_dir, cancer_type, "/cnv/brca_cnvw_cancer.Rdata", sep=""))
load(paste(parsed_data_dir, cancer_type, "/expr/brca_expr_cancer.Rdata", sep=""))
load(paste(parsed_data_dir, cancer_type, "/muta/brca_muta_cancer.Rdata", sep=""))

# TODO: Add other data sets

# Subset data with nice_genes --------------------------------------------------
cat("Subsetting data..\n")
brca_cnvw_cancer_min <- brca_cnvw_cancer[genename %in% nice_genes]
brca_expr_cancer_min <- brca_expr_cancer[gene %in% nice_genes]
brca_muta_cancer_min <- brca_muta_cancer

# Save to database -------------------------------------------------------------

cat("Saving to database..\n")
db <- dbConnect(SQLite(), dbname=paste(parsed_data_dir, cancer_type, "/",
                                       cancer_type, "_min.sqlite", sep=""))

if(length(dbListTables(db)) != 0){
  stop("Database isn't empty\n")
}

dbWriteTable(conn = db, name = "nice_genes",
             value = data.frame(nice_genes), row.names=F)
dbWriteTable(conn = db, name = "brca_cnvw_cancer",
             value = brca_cnvw_cancer_min, row.names = F)
dbWriteTable(conn = db, name = "brca_expr_cancer",
             value = brca_expr_cancer_min, row.names = F)
dbWriteTable(conn = db, name = "brca_muta_cancer",
             value = brca_muta_cancer_min, row.names = F)

dbDisconnect(db)
