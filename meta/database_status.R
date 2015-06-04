library(RSQLite)

raw_data_dir    <- "/share/scratch/arj32/raw_data/"
parsed_data_dir <- "/share/scratch/arj32/parsed_data/"

cancer_types <- c("brca", "ov", "paad", "prad")
data_types   <- c("clin", "cnv", "expr", "meth", "muta", "mirn")

status_list <- list()

for (c_type in cancer_types) {
  raw   <- rep(" ", length(data_types))
  Rdata <- rep(" ", length(data_types))
  sql   <- rep(" ", length(data_types))
  names(raw)   <- data_types
  names(Rdata) <- data_types
  names(sql)   <- data_types

  for (d_type in data_types) {
    # Check files in raw_data --------------------------------------------------
    if (d_type == "clin"){
      check_file <- paste(raw_data_dir, c_type,
                          "/clin/nationwidechildrens.org_clinical_patient_",
                          c_type, ".txt", sep="")
      if (file.exists(check_file)) {
        raw[d_type] <- "x"
      }
    } else {
      check_file <- paste(raw_data_dir, c_type, "/", d_type,
                          "/FILE_SAMPLE_MAP.txt", sep="")
      if (file.exists(check_file)) {
        raw[d_type] <- "x"
      }
    }

    # Check Rdata files in parsed_data -----------------------------------------
    if (d_type == "cnv") {
      check_file_c <- paste(parsed_data_dir, c_type, "/", d_type, "/", c_type,
                            "_", "cnvw_cancer.Rdata", sep="")
      check_file_n <- paste(parsed_data_dir, c_type, "/", d_type, "/", c_type,
                            "_", "cnvw_normal.Rdata", sep="")
    } else {
      check_file_c <- paste(parsed_data_dir, c_type, "/", d_type, "/", c_type,
                            "_", d_type, "_cancer.Rdata", sep="")
      check_file_n <- paste(parsed_data_dir, c_type, "/", d_type, "/", c_type,
                            "_", d_type, "_normal.Rdata", sep="")
    }
    if (file.exists(check_file_c) && file.exists(check_file_n)) {
      Rdata[d_type] <- "b"
    } else if (file.exists(check_file_c) && !file.exists(check_file_n)) {
      Rdata[d_type] <- "c"
    } else if (!file.exists(check_file_c) && file.exists(check_file_n)) {
      Rdata[d_type] <- "n"
    }

    # Check SQLite database in parsed_data -------------------------------------
    sql_file <- paste(parsed_data_dir, c_type, "/", c_type, ".sqlite", sep="")
    if (file.exists(sql_file)){
      db <- dbConnect(SQLite(), dbname=sql_file)
      tables <- dbListTables(db)
      dbDisconnect(db)
      for (d_type in data_types) {
        if (d_type == "cnv") {
          check_c <- any(grepl(paste("cnvw_cancer", sep=""), tables))
          check_n <- any(grepl(paste("cnvw_normal", sep=""), tables))
        } else {
          check_c <- any(grepl(paste(d_type, "_cancer", sep=""), tables))
          check_n <- any(grepl(paste(d_type, "_normal", sep=""), tables))
        }
        if (check_c && check_n) {
          sql[d_type] <- "b"
        } else if (check_c && !check_n) {
          sql[d_type] <- "c"
        } else if (!check_c && check_n) {
          sql[d_type] <- "n"
        }
      }
    }
  }
  status_list[[c_type]] <- data_frame(raw   = paste("[ ", raw,   " ]", sep=""),
                                      Rdata = paste("[ ", Rdata, " ]", sep=""),
                                      sql   = paste("[ ", sql,   " ]", sep=""),
                                      row.names = data_types)
}
# Output status to console -----------------------------------------------------
print(status_list)
