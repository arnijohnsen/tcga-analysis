library(data.table)
library(RSQLite)

cancer_type <- "brca"
raw_data_dir    <- "/share/scratch/arj32/raw_data/"
parsed_data_dir <- "/share/scratch/arj32/parsed_data/"

# Open connection to database --------------------------------------------------
db <- dbConnect(SQLite(), dbname=paste(parsed_data_dir, cancer_type, "/",
                                       cancer_type, ".sqlite", sep=""))

# Check if db contains anything ------------------------------------------------
if(length(dbListTables(db)) != 0){
  stop("Database isn't empty\n")
}

# Read data.tables one by one and add them to SQLite database ------------------
data_types <- c("cnvw", "expr", "meth", "muta", "mirn")
data_dirs  <- c("cnv",  "expr", "meth", "muta", "mirn")
for (i in 1:length(data_types)){
  cat("Adding", data_types[i], "data..\n")
  for (j in c("cancer", "normal")){
    if(!(data_types[i] == "muta" & j == "normal")){
      cat("..", j, "\n", sep="")
      cat("....reading\n")
      load(paste(parsed_data_dir, cancer_type, "/", data_dirs[i], "/",
                 cancer_type, "_", data_types[i], "_", j, ".Rdata", sep=""))
      cat("....writing\n")
      dbWriteTable(conn  = db,
                   name  = paste(data_types[i], "_", j, sep=""),
                   value = get(paste(cancer_type, data_types[i], j, sep=".")),
                   row_names = F)
      rm(list=paste(cancer_type, data_types[i], j, sep="."))
      gc()
    }
  }
}

# Close connection to database -------------------------------------------------
dbDisconnect(db)
