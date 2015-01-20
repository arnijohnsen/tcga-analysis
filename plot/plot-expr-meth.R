# TODO:
#  x Enable plotting of body and enchancers
#  x Read search string from command line
#  x Enable search for probes
#  x Don't reload files
#  - Allow user to select which types (enhancer, promoter,..) are used
#  x Load only needed files
#  x Do't plot empty plots
#  - Implement new expression reading scripts and datafiles (should be simple)
#  - Enable abline
#  - Add regression statistics to plot
#  - Add annotation to plot
#  - Enable option to show only matches samples and draw line between matching
#  - Use ggplot2 instead of plot
#  - Write this as a function maybe?

# Should be load() with Rdata file
if(!exists("linked.probes.genes")){
  linked.probes.genes <- read.table("../Rdata/BRCA-new/probe-annotation/linked-probes-genes.txt", stringsAsFactors=F)
}

cat("Loading expression\n")
if(!exists("BRCA.NEA")){
  load("../Rdata/BRCA/data/BRCA-NEA.Rdata")
}
if(!exists("BRCA.CEA")){
  load("../Rdata/BRCA/data/BRCA-CEA.Rdata")
}

cat("Loading methylation\n")
for(set in c("mn","mc","uc")){
  for(type in c("bi", "bs", "bn", "en", "pr")){
    if(!exists(paste(set,".",type,sep=""))){
      load(paste("../parsed-data/BRCA/methylation/", set, "-", type, ".Rdata", sep=""))
    }
  }
}

if(!exists("mn.list")){
  mn.list <- list(body.island=mn.bi, body.shore=mn.bs, body.none=mn.bn, enhancer=mn.en, promoter=mn.pr)
  mc.list <- list(body.island=mc.bi, body.shore=mc.bs, body.none=mc.bn, enhancer=mc.en, promoter=mc.pr)
  uc.list <- list(body.island=uc.bi, body.shore=uc.bs, body.none=uc.bn, enhancer=uc.en, promoter=uc.pr)
}

search <- readline("Enter gene or probe name: ")
search <- paste("^", search, "$", sep="")

idx <- c(grep(search, linked.probes.genes$gene), grep(search, linked.probes.genes$probe))

for(i in idx){
  probe  <- linked.probes.genes$probe[i]
  gene   <- linked.probes.genes$gene[i]
  status <- linked.probes.genes$status[i]
  cat("Probe:", probe, "\tgene:", gene, "\tstatus: ", status, "\n")

  # Matched normal
  x.mn.all <- mn.list[[status]][[probe]]
  names(x.mn.all) <- substring(rownames(mn.list[[status]]),1,14)
  y.mn.all <- BRCA.NEA[[gene]]
  names(y.mn.all) <- rownames(BRCA.NEA)
  mn.samples <- intersect(names(x.mn.all),names(y.mn.all))

  x.mn.plot <- x.mn.all[mn.samples]
  y.mn.plot <- y.mn.all[mn.samples]

  # Matched cancer
  x.mc.all <- mc.list[[status]][[probe]]
  names(x.mc.all) <- substring(rownames(mc.list[[status]]),1,14)
  y.mc.all <- BRCA.CEA[[gene]]
  names(y.mc.all) <- rownames(BRCA.CEA)
  mc.samples <- intersect(names(x.mc.all),names(y.mc.all))

  x.mc.plot <- x.mc.all[mc.samples]
  y.mc.plot <- y.mc.all[mc.samples]

  # Unmatched cancer
  x.uc.all <- uc.list[[status]][[probe]]
  names(x.uc.all) <- substring(rownames(uc.list[[status]]),1,14)
  y.uc.all <- BRCA.CEA[[gene]]
  names(y.uc.all) <- rownames(BRCA.CEA)
  uc.samples <- intersect(names(x.uc.all),names(y.uc.all))

  x.uc.plot <- x.uc.all[uc.samples]
  y.uc.plot <- y.uc.all[uc.samples]

  ymax <- max(c(y.mn.plot, y.mc.plot, y.uc.plot))

  if(!all(is.na(c(x.mn.plot,x.mc.plot,x.uc.plot)))){
    plot(x.mn.plot, y.mn.plot, xlim=c(0,100), ylim=c(0,ymax), col=rgb(0,0,1,0.5), pch=20)
    points(x.mc.plot, y.mc.plot, col=rgb(0,1,0,0.5), pch=20)
    points(x.uc.plot, y.uc.plot, col=rgb(1,0,0,0.5), pch=20)

    cat("Press [enter] for next plot")
    line <- readline()
  } else {
    cat("All methylation values are NA\n")
  }
}
