# TODO:
#  - Enable plotting of body and enchancers
#  - Read search string from command line
#  - Enable search for probes
#  - Don't reload files
#  - Load only needed files
#  - Implement new expression reading script (should be simple)
#  - Enable abline
#  - Add regression statistics to plot
#  - Add annotation to plot
#  - Enable option to show only matches samples and draw line between matching
#  - Use ggplot2 instead of plot
#  - Write this as a function maybe?

# Should be load() with Rdata file
linked.probes.genes <- read.table("../Rdata/BRCA-new/probe-annotation/linked-probes-genes.txt", stringsAsFactors=F)

# Should be read from command line
search <- "BRCA" 

# Should support probe search
idx <- grep(search, linked.probes.genes$gene)

cat("Loading expression\n")
load("../Rdata/BRCA/data/BRCA-NEA.Rdata")
load("../Rdata/BRCA/data/BRCA-CEA.Rdata")

cat("Loading methylation\n")
for(set in c("mn","mc","uc")){
  for(type in c("bi", "bs", "bn", "en", "pr")){
    load(paste("../parsed-data/BRCA/methylation/", set, "-", type, ".Rdata", sep=""))
  }
}

for(i in idx){
  probe  <- linked.probes.genes$probe[i]
  gene   <- linked.probes.genes$gene[i]
  status <- linked.probes.genes$status[i]
  cat("Probe:", probe, "\tgene:", gene, "\tstatus: ", status, "\n")

  if(status=="promoter"){
    # Matched normal
    x.mn.all <- mn.pr[[probe]]
    names(x.mn.all) <- substring(rownames(mn.pr),1,14)
    y.mn.all <- BRCA.NEA[[gene]]
    names(y.mn.all) <- rownames(BRCA.NEA)
    mn.samples <- intersect(names(x.mn.all),names(y.mn.all))

    x.mn.plot <- x.mn.all[mn.samples]
    y.mn.plot <- y.mn.all[mn.samples]

    # Matched cancer
    x.mc.all <- mc.pr[[probe]]
    names(x.mc.all) <- substring(rownames(mc.pr),1,14)
    y.mc.all <- BRCA.CEA[[gene]]
    names(y.mc.all) <- rownames(BRCA.CEA)
    mc.samples <- intersect(names(x.mc.all),names(y.mc.all))

    x.mc.plot <- x.mc.all[mc.samples]
    y.mc.plot <- y.mc.all[mc.samples]

    # Unmatched cancer
    x.uc.all <- uc.pr[[probe]]
    names(x.uc.all) <- substring(rownames(uc.pr),1,14)
    y.uc.all <- BRCA.CEA[[gene]]
    names(y.uc.all) <- rownames(BRCA.CEA)
    uc.samples <- intersect(names(x.uc.all),names(y.uc.all))

    x.uc.plot <- x.uc.all[uc.samples]
    y.uc.plot <- y.uc.all[uc.samples]

    ymax <- max(c(y.mn.plot, y.mc.plot, y.uc.plot))

    plot(x.mn.plot, y.mn.plot, xlim=c(0,100), ylim=c(0,ymax), col=rgb(0,0,1,0.5), pch=20)
    points(x.mc.plot, y.mc.plot, col=rgb(0,1,0,0.5), pch=20)
    points(x.uc.plot, y.uc.plot, col=rgb(1,0,0,0.5), pch=20)

    cat("Press [enter] for next plot")
    line <- readline()
  }
}
