data.file.dir = "../rawdata/BRCA/Methylation/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"
output.dir    = "../parsed-data/BRCA/methylation/"

load("../parsed-data/BRCA/annotation/methylation-files.Rdata")
load("../parsed-data/BRCA/annotation/sorted-probes.Rdata")

##########################
# Matched normal samples #
##########################

write(paste(sorted.probes$bi, collapse="\t"), file=paste(output.dir, "mn-bi.tmp", sep=""))
write(paste(sorted.probes$bs, collapse="\t"), file=paste(output.dir, "mn-bs.tmp", sep=""))
write(paste(sorted.probes$bn, collapse="\t"), file=paste(output.dir, "mn-bn.tmp", sep=""))
write(paste(sorted.probes$en, collapse="\t"), file=paste(output.dir, "mn-en.tmp", sep=""))
write(paste(sorted.probes$pr, collapse="\t"), file=paste(output.dir, "mn-pr.tmp", sep=""))

n.mn <- length(mn.files)
for(i in 1:n.mn){
  cat("Reading matched normal file", i, "of", n.mn, "\n")
  raw.file.data <- read.table(paste(data.file.dir, mn.files[i], sep=""), header=TRUE, sep="\t", quote="\"", skip=1, stringsAsFactors=F)
  sample.name <- sub(".txt","", sub(".*TCGA","TCGA", mn.files[i]))

  tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% sorted.probes$bi]
  write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir, "mn-bi.tmp", sep=""), append=T)

  tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% sorted.probes$bs]
  write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir, "mn-bs.tmp", sep=""), append=T)

  tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% sorted.probes$bn]
  write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir, "mn-bn.tmp", sep=""), append=T)

  tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% sorted.probes$en]
  write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir, "mn-en.tmp", sep=""), append=T)

  tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% sorted.probes$pr]
  write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir, "mn-pr.tmp", sep=""), append=T)
}

##########################
# Matched cancer samples #
##########################

write(paste(sorted.probes$bi, collapse="\t"), file=paste(output.dir, "mc-bi.tmp", sep=""))
write(paste(sorted.probes$bs, collapse="\t"), file=paste(output.dir, "mc-bs.tmp", sep=""))
write(paste(sorted.probes$bn, collapse="\t"), file=paste(output.dir, "mc-bn.tmp", sep=""))
write(paste(sorted.probes$en, collapse="\t"), file=paste(output.dir, "mc-en.tmp", sep=""))
write(paste(sorted.probes$pr, collapse="\t"), file=paste(output.dir, "mc-pr.tmp", sep=""))

n.mc <- length(mc.files)
for(i in 1:n.mc){
  cat("Reading matched normal file", i, "of", n.mc, "\n")
  raw.file.data <- read.table(paste(data.file.dir, mc.files[i], sep=""), header=TRUE, sep="\t", quote="\"", skip=1, stringsAsFactors=F)
  sample.name <- sub(".txt","", sub(".*TCGA","TCGA", mc.files[i]))

  tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% sorted.probes$bi]
  write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir, "mc-bi.tmp", sep=""), append=T)

  tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% sorted.probes$bs]
  write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir, "mc-bs.tmp", sep=""), append=T)

  tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% sorted.probes$bn]
  write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir, "mc-bn.tmp", sep=""), append=T)

  tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% sorted.probes$en]
  write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir, "mc-en.tmp", sep=""), append=T)

  tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% sorted.probes$pr]
  write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir, "mc-pr.tmp", sep=""), append=T)
}

############################
# Unmatched cancer samples #
############################

write(paste(sorted.probes$bi, collapse="\t"), file=paste(output.dir, "uc-bi.tmp", sep=""))
write(paste(sorted.probes$bs, collapse="\t"), file=paste(output.dir, "uc-bs.tmp", sep=""))
write(paste(sorted.probes$bn, collapse="\t"), file=paste(output.dir, "uc-bn.tmp", sep=""))
write(paste(sorted.probes$en, collapse="\t"), file=paste(output.dir, "uc-en.tmp", sep=""))
write(paste(sorted.probes$pr, collapse="\t"), file=paste(output.dir, "uc-pr.tmp", sep=""))

n.uc <- length(uc.files)
for(i in 1:n.uc){
  cat("Reading matched normal file", i, "of", n.uc, "\n")
  raw.file.data <- read.table(paste(data.file.dir, uc.files[i], sep=""), header=TRUE, sep="\t", quote="\"", skip=1, stringsAsFactors=F)
  sample.name <- sub(".txt","", sub(".*TCGA","TCGA", uc.files[i]))

  tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% sorted.probes$bi]
  write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir, "uc-bi.tmp", sep=""), append=T)

  tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% sorted.probes$bs]
  write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir, "uc-bs.tmp", sep=""), append=T)

  tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% sorted.probes$bn]
  write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir, "uc-bn.tmp", sep=""), append=T)

  tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% sorted.probes$en]
  write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir, "uc-en.tmp", sep=""), append=T)

  tmp <- raw.file.data$Beta_value[raw.file.data$Composite.Element.REF %in% sorted.probes$pr]
  write(paste(sample.name, paste(round(tmp*100,digits=0),collapse="\t"),sep="\t"), file=paste(output.dir, "uc-pr.tmp", sep=""), append=T)
}
quit(save="no")
