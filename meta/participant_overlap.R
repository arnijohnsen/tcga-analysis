library(data.table)
library(limma)

# Define cancer type, raw and parsed data directories --------------------------
cancer_type <- "brca"
raw_data_dir    <- "/share/scratch/arj32/raw_data/"
parsed_data_dir <- "/share/scratch/arj32/parsed_data/"

# Read pariticipant list for cnv, expr and meth. TODO: add muta ----------------
cnv <- fread(paste(parsed_data_dir, cancer_type,
                           "/info/cnv_participants_txt", sep=""))
expr <- fread(paste(parsed_data_dir, cancer_type,
                           "/info/expr_participants_txt", sep=""))
meth <- fread(paste(parsed_data_dir, cancer_type,
                           "/info/meth_participants_txt", sep=""))
muta <- fread(paste(parsed_data_dir, cancer_type,
                           "/info/muta_participants_txt", sep=""))
mirn <- fread(paste(parsed_data_dir, cancer_type,
                           "/info/mirn_participants_txt", sep=""))

cnv_c  <-  cnv[!is_na(cancer_barcode)]
expr_c <- expr[!is_na(cancer_barcode)]
meth_c <- meth[!is_na(cancer_barcode)]
muta_c <- muta[!is_na(cancer_barcode)]
mirn_c <- mirn[!is_na(cancer_barcode)]

cnv_n  <-  cnv[!is_na(normal_barcode)]
expr_n <- expr[!is_na(normal_barcode)]
meth_n <- meth[!is_na(normal_barcode)]
muta_n <- muta[!is_na(normal_barcode)]
mirn_n <- mirn[!is_na(normal_barcode)]

cnv_b  <-  cnv[!is_na(cancer_barcode)&!is_na(normal_barcode)]
expr_b <- expr[!is_na(cancer_barcode)&!is_na(normal_barcode)]
meth_b <- meth[!is_na(cancer_barcode)&!is_na(normal_barcode)]
muta_b <- muta[!is_na(cancer_barcode)&!is_na(normal_barcode)]
mirn_b <- mirn[!is_na(cancer_barcode)&!is_na(normal_barcode)]

all_part <- unique(c(cnv$participant,
                     expr$participant,
                     meth$participant,
                     muta$participant,
                     mirn$participant))

DF.c <- data.table(cnv  = all_part %in% cnv_c$participant,
                   expr = all_part %in% expr_c$participant,
                   meth = all_part %in% meth_c$participant,
                   muta = all_part %in% muta_c$participant,
                   mirn = all_part %in% mirn_c$participant)

DF.n <- data.table(cnv  = all_part %in% cnv_n$participant,
                   expr = all_part %in% expr_n$participant,
                   meth = all_part %in% meth_n$participant,
                   muta = all_part %in% muta_n$participant,
                   mirn = all_part %in% mirn_n$participant)

DF.b <- data.table(cnv  = all_part %in% cnv_b$participant,
                   expr = all_part %in% expr_b$participant,
                   meth = all_part %in% meth_b$participant,
                   muta = all_part %in% muta_b$participant,
                   mirn = all_part %in% mirn_b$participant)

VC.c <- vennCounts(DF.c)
VC.n <- vennCounts(DF.n)
VC.b <- vennCounts(DF.b)
overlap <- data.table(cbind(VC.c, VC.n[,6], VC.b[,6]))
setnames(overlap, c("cnv", "expr", "meth", "muta", "mirn", "c", "n", "b"))

n <- dim(overlap)[1]

overlap$c_sum <- rep(0, n)
overlap$n_sum <- rep(0, n)
overlap$b_sum <- rep(0, n)

for (i in 1:n){
  idx <- sapply(1:32, function(x){
    return(all(overlap[x,1:5,with=F] >= overlap[i,1:5,with=F]))
  })
  overlap$c_sum[i] <- sum(overlap$c[idx])
  overlap$n_sum[i] <- sum(overlap$n[idx])
  overlap$b_sum[i] <- sum(overlap$b[idx])
}

#sink(file=paste(parsed_data_dir, cancer_type,
#                "/info/participant_overlap.txt", sep=""))
#print(overlap)
#sink()
