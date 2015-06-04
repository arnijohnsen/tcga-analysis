library(data.table) # Uses fread (much faster than read.table)

load("../parsed-data/BRCA/annotation/methylation-files.Rdata")
load("../parsed-data/BRCA/annotation/sorted-probes.Rdata")

barcode <- function(x){
  # Extracts barcode from file name (.txt files)
  substr(x, nchar(x)-31, nchar(x)-4)
}

for(set in c("mn","mc","uc")){
  for(type in c("bi", "bs", "bn", "en", "pr")){
    cat("Reading",paste("../parsed-data/BRCA/methylation/", set, "-", type, ".tmp", sep=""),"..\n")

    tmp <- fread(paste("../parsed-data/BRCA/methylation/", set, "-", type, ".tmp", sep=""), stringsAsFactors=F, data.table=F)
    tmp$V1 <- NULL
    colnames(tmp) <- sorted.probes[[type]]
    rownames(tmp) <- barcode(get(paste(set,".files",sep="")))

    df.name <- paste(set,".",type,sep="")
    assign(df.name, tmp)

    cat("Saving", df.name, "as", paste("../parsed-data/BRCA/methylation/", set, "-", type, ".Rdata", sep=""), "..\n")

    save(list=df.name,  file=paste("../parsed-data/BRCA/methylation/", set, "-", type, ".Rdata", sep=""))
    rm(list=df.name)
  }
}
quit(save="no")
