library(RSQLite)

raw.data.dir    <- "/share/scratch/arj32/raw-data/"
parsed.data.dir <- "/share/scratch/arj32/parsed-data/"

cancer.types <- c("brca", "ov", "paad", "prad")
data.types   <- c("clin", "cnv", "expr", "meth", "muta", "mirn")

status.list <- list()

for (c.type in cancer.types) {
  raw    <- rep(FALSE, length(data.types))
  Rdata  <- rep(FALSE, length(data.types))
  sqlite <- rep(FALSE, length(data.types))
  names(raw)    <- data.types
  names(Rdata)  <- data.types
  names(sqlite) <- data.types

  for (d.type in data.types) {
    # Check files in raw-data --------------------------------------------------
    raw[d.type] <- file.exists(paste(raw.data.dir, c.type, "/", d.type,
                                 "/README_DCC.txt", sep=""))

    # Check Rdata files in parsed-data -----------------------------------------
    Rdata[d.type] <- file.exists(paste(parsed.data.dir, c.type, "/", d.type,
                                       "/", c.type, "-", d.type,
                                       "-cancer.Rdata", sep=""))
    if (d.type == "cnv") {
      Rdata[d.type] <- file.exists(paste(parsed.data.dir, c.type, "/", d.type,
                                         "/", c.type, "-", "cnvw",
                                         "-cancer.Rdata", sep=""))
    }
  }
  # Check SQLite database in parsed-data -------------------------------------
  sql.file <- paste(parsed.data.dir, c.type, "/", c.type, ".sqlite", sep="")
  if (file.exists(sql.file)){
    db <- dbConnect(SQLite(), dbname=sql.file)
    tables <- dbListTables(db)
    dbDisconnect(db)
    for (d.type in data.types) {
      sqlite[d.type] <- any(grepl(d.type, tables))
    }
  }
  status.list[[c.type]] <- data.frame(raw, Rdata, sqlite)
}
print(status.list)
