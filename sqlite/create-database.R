library(data.table)
library(RSQLite)

cancer.type <- "brca"
raw.data.dir    <- "/share/scratch/arj32/raw-data/"
parsed.data.dir <- "/share/scratch/arj32/parsed-data/"

# Open database
db <- dbConnect(SQLite(), dbname=paste(parsed.data.dir, cancer.type, "/", cancer.type, ".sqlite", sep=""))

# Check if db contains anything
if(length(dbListTables(db)) != 0){
  stop("Database isn't empty\n")
}

data.types <- c("cnvw", "expr", "meth")
data.dirs  <- c("cnv",  "expr", "meth")
for(i in 1:3){
  cat("Adding", data.types[i], "data..\n")
  cat("..normal\n")
  cat("....reading\n")
  load(paste(parsed.data.dir, cancer.type, "/", data.dirs[i], "/", cancer.type, "-", data.types[i], "-normal.Rdata", sep=""))
  cat("....writing\n")
  dbWriteTable(conn  = db,
               name  = paste(data.types[i], "_normal", sep=""),
               value = get(paste(cancer.type, data.types[i], "normal", sep=".")),
               row.names = F)
  rm(list=paste(cancer.type, data.types[i], "normal", sep="."))
}

#cat(dbListTables(db))
dbDisconnect(db)
