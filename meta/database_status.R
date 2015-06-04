library(RSQLite)

raw.data.dir    <- "/share/scratch/arj32/raw_data/"
parsed.data.dir <- "/share/scratch/arj32/parsed_data/"

cancer.types <- c("brca", "ov", "paad", "prad")
data.types   <- c("clin", "cnv", "expr", "meth", "muta", "mirn")

status.list <- list()

for (c.type in cancer.types) {
  raw   <- rep(" ", length(data.types))
  Rdata <- rep(" ", length(data.types))
  sql   <- rep(" ", length(data.types))
  names(raw)   <- data.types
  names(Rdata) <- data.types
  names(sql)   <- data.types

  for (d.type in data.types) {
    # Check files in raw_data --------------------------------------------------
    if (d.type == "clin"){
      check.file <- paste(raw.data.dir, c.type,
                          "/clin/nationwidechildrens.org_clinical_patient_",
                          c.type, ".txt", sep="")
      if (file.exists(check.file)) {
        raw[d.type] <- "x"
      }
    } else {
      check.file <- paste(raw.data.dir, c.type, "/", d.type,
                          "/FILE_SAMPLE_MAP.txt", sep="")
      if (file.exists(check.file)) {
        raw[d.type] <- "x"
      }
    }

    # Check Rdata files in parsed_data -----------------------------------------
    if (d.type == "cnv") {
      check.file.c <- paste(parsed.data.dir, c.type, "/", d.type, "/", c.type,
                            "_", "cnvw_cancer.Rdata", sep="")
      check.file.n <- paste(parsed.data.dir, c.type, "/", d.type, "/", c.type,
                            "_", "cnvw_normal.Rdata", sep="")
    } else {
      check.file.c <- paste(parsed.data.dir, c.type, "/", d.type, "/", c.type,
                            "_", d.type, "_cancer.Rdata", sep="")
      check.file.n <- paste(parsed.data.dir, c.type, "/", d.type, "/", c.type,
                            "_", d.type, "_normal.Rdata", sep="")
    }
    if (file.exists(check.file.c) && file.exists(check.file.n)) {
      Rdata[d.type] <- "b"
    } else if (file.exists(check.file.c) && !file.exists(check.file.n)) {
      Rdata[d.type] <- "c"
    } else if (!file.exists(check.file.c) && file.exists(check.file.n)) {
      Rdata[d.type] <- "n"
    }

    # Check SQLite database in parsed_data -------------------------------------
    sql.file <- paste(parsed.data.dir, c.type, "/", c.type, ".sqlite", sep="")
    if (file.exists(sql.file)){
      db <- dbConnect(SQLite(), dbname=sql.file)
      tables <- dbListTables(db)
      dbDisconnect(db)
      for (d.type in data.types) {
        if (d.type == "cnv") {
          check.c <- any(grepl(paste("cnvw_cancer", sep=""), tables))
          check.n <- any(grepl(paste("cnvw_normal", sep=""), tables))
        } else {
          check.c <- any(grepl(paste(d.type, "_cancer", sep=""), tables))
          check.n <- any(grepl(paste(d.type, "_normal", sep=""), tables))
        }
        if (check.c && check.n) {
          sql[d.type] <- "b"
        } else if (check.c && !check.n) {
          sql[d.type] <- "c"
        } else if (!check.c && check.n) {
          sql[d.type] <- "n"
        }
      }
    }
  }
  status.list[[c.type]] <- data.frame(raw   = paste("[ ", raw,   " ]", sep=""),
                                      Rdata = paste("[ ", Rdata, " ]", sep=""),
                                      sql   = paste("[ ", sql,   " ]", sep=""),
                                      row.names = data.types)
}
# Output status to console -----------------------------------------------------
print(status.list)
