library(data.table)
library(RSQLite)

cancer.type <- "brca"
raw.data.dir    <- "/share/scratch/arj32/raw-data/"
parsed.data.dir <- "/share/scratch/arj32/parsed-data/"

# Open connection to database --------------------------------------------------
db <- dbConnect(SQLite(), dbname=paste(parsed.data.dir, cancer.type, "/",
                                       cancer.type, ".sqlite", sep=""))

# Check if db contains anything ------------------------------------------------
if(length(dbListTables(db)) != 0){
  stop("Database isn't empty\n")
}

# Read data.tables one by one and add them to SQLite database ------------------
data.types <- c("cnvw", "expr", "meth", "muta", "mirn")
data.dirs  <- c("cnv",  "expr", "meth", "muta", "mirn")
for (i in 1:length(data.types)){
  cat("Adding", data.types[i], "data..\n")
  for (j in c("cancer", "normal")){
    if(!(data.types[i] == "muta" & j == "normal")){
      cat("..", j, "\n", sep="")
      cat("....reading\n")
      load(paste(parsed.data.dir, cancer.type, "/", data.dirs[i], "/",
                 cancer.type, "-", data.types[i], "-", j, ".Rdata", sep=""))
      cat("....writing\n")
      dbWriteTable(conn  = db,
                   name  = paste(data.types[i], "_", j, sep=""),
                   value = get(paste(cancer.type, data.types[i], j, sep=".")),
                   row.names = F)
      rm(list=paste(cancer.type, data.types[i], j, sep="."))
      gc()
    }
  }
}

# Close connection to database -------------------------------------------------
dbDisconnect(db)
