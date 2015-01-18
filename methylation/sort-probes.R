cat("Reading probes annotation file..\n")
probe.annotations <- read.table("../rawdata/AnnotationFiles/GenomeStudioProbeAnnotations.txt", header=TRUE, sep="\t", quote="\"", stringsAsFactors=FALSE)

cat("Setting probe status..\n")
probe.annotations$status <- 0

cat("  ..gene bodies\n")
is.body   <- !grepl("TSS200|TSS1500|1stExon|5'UTR|3'UTR|^$", probe.annotations$UCSC_REFGENE_GROUP)
is.island <-  grepl("Island", probe.annotations$RELATION_TO_UCSC_CPG_ISLAND)
is.shore  <-  grepl("Shore",  probe.annotations$RELATION_TO_UCSC_CPG_ISLAND)
is.none   <- !is.island & !is.shore
probe.annotations$status[is.body & is.island] <- "body.island"
probe.annotations$status[is.body & is.shore]  <- "body.shore"
probe.annotations$status[is.body & is.none]   <- "body.none"

cat("  ..enhancers\n")
probe.annotations$status[probe.annotations$ENHANCER] <- "enhancer"
cat("  ..promoters\n")
probe.annotations$status[grepl("TSS200|5'UTR", probe.annotations$UCSC_REFGENE_GROUP)] <- "promoter"



body.island.probes <- probe.annotations$TargetID[probe.annotations$status == "body.island"]
body.shore.probes  <- probe.annotations$TargetID[probe.annotations$status == "body.shore" ]
body.none.probes   <- probe.annotations$TargetID[probe.annotations$status == "body.none"  ]
enhancer.probes    <- probe.annotations$TargetID[probe.annotations$status == "enhancer"   ]
promoter.probes    <- probe.annotations$TargetID[probe.annotations$status == "promoter"   ]

sorted.probes <- list(bi = body.island.probes, bs = body.shore.probes, bn = body.none.probes, en = enhancer.probes, pr = promoter.probes)

# Save as .Rdata
save(sorted.probes, file="../parsed-data/BRCA/methylation/sorted-probes.txt")
quit(save="no")
