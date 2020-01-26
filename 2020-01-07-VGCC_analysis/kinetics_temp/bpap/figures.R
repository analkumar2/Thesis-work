source("figures-common.R")

dataset.cache <- TRUE

## Set colours of lines
palette(c("black","gray","red","blue","green","orange"))

colwidth <- 3.27                        # Width of column

## Function to plot everything
plot.all.features.dist.size <- function(group="car-mag-new") {
  panlabel.adj <- -5

  ######################################################################
  ## Figure 1
  ##
  ## Calcium concentration during single BAP as a function of
  ## distance along branch
  ######################################################################
  
  dataset <-  "ca1_poirazi-somainj-Ra050-cv1"
  r <- get.dataset(dataset, dir=file.path("datastore", group))

  figmessage("Figure 1: Calcium traces")
  ## Panel C: Traces of V, Ca in trunk and Ca in spines
  stdpostscript(file=paste(dataset, "-synrecs", sep=""), group=group,
                width=0.6*colwidth, height=4.3)
  # plot.synrecs(r, panlabs=c("C", "", "", ""),  xlim.closeup=NULL)
  par(mfrow=c(4,1))
  stdpars()
  par(oma=c(0,0,0.5,0))
  plot.synrec.closeup(r$tsoma, r$screc_vsyn,
                      xlab=expression(paste(italic(t), " (ms)")),
                      ylab=expression(paste(italic(V), " (mV)")),
                      panlab="",
                      x.extra=r$tsoma, y.extra=r$vsoma,
                      ylim=c(-80,20),
                      xlim.closeup=NULL)
  panlabel("E", line=0.2)
  
  plot.synrec.closeup(r$tsoma,1000*r$screc_catree,
                      xlab=expression(paste(italic(t), " (ms)")),
                      ylab=expression(paste("[Ca] (", mu,"M)",sep="")),
                      panlab="",
                      xlim.closeup=NULL,
                      ylim=c(0,2))
  panlabel("F", line=0.2)
  ## plot.synrec.closeup(r$tsoma,1000*r$screc_casyn,
  ##                     xlab=expression(paste(italic(t), " (ms)")),
  ##                     ylab=expression(paste("Syn. [Ca] (", mu,"M)",sep="")),
  ##                     panlab="",
  ##                     col=c(3,4,5,6),
  ##                     lwd=2,lty=c(1,1,1),
  ##                     xlim.closeup=NULL)
  ## dev.off()
  
  ## stdpostscript(file=dataset, group=group, width=2.5, height=3.3)
  ## par(mfrow=c(2, 1))
  par(mar=stdparlist$mar)

  figmessage("Figure 1D: Membrane potential")
  ## Bottom row: Membrane potential
  vrest <- -70
  plot.feature(r, x=r$distances,
               feature=r$vtreemax_mean-vrest,
               xlab=dist.lab, xlim=dist.lim,
               ylim=c(0, 100),
               fit=FALSE,
               ylab=expression("Peak m.p. (mV)"),
               legend.pos="topleft",
               plotinds=1:length(r$vtreemax_mean),
               yaxp=c(0, 100, 2), ordinate.dependent=FALSE)
  panlabel("G", line=0.2)

  figmessage("Figure 1E: Calcium in trunk")
  plot.feature(r, x=r$distances,
               feature=r$catreemax_mean*1000,
               xlab=dist.lab, xlim=dist.lim,
               ylim=c(0, 0.002*1000),
               fit=FALSE,
               ylab=expression(paste("Peak [Ca]", " (", mu, "M)",sep="")),
               legend.pos="topleft",
               #plotinds=1:length(r$catreemax_mean),
               plotinds=r$trunk_sri+1, ordinate.dependent=FALSE)
  panlabel("H", line=0.2)
  dev.off()

  ######################################################################
  ##
  ## Figure 2
  ##
  ######################################################################
  
  dataset <-  "ca1_poirazi-dendburst-s240-j00t1-a200-n45-bv-r170-sc0-Ra050-nr0100-cv1"
  r <- get.dataset(dataset, dir=file.path("datastore", group))

  ## Save for ks test later
  r.sc0 <- r
  
  ## Figure 2 - time courses from spines
  figmessage("Figure 2B-F: Time courses in spines")
  stdpostscript(file=paste(dataset, "-synrecs", sep=""), group=group,
                width=colwidth, height=4)
  par(mar=c(1.5, 2.5, 0.1, 0.8))
  plot.synrecs(r, panlabs=c("B", "C", "D", "E", "F"), xlim.closeup=c(22, 28),
               show.IR=TRUE)
  dev.off()

  ## Figure 2 - Measures of distance with SC stimulation
  figmessage("Figure 2G-L: Measures of distance with SC stimulation")
  stdpostscript(file=dataset, group=group, width=1.5*colwidth, height=3)
  par(mfrow=c(2, 3))
  stdpars()
  
  ## Top row: peak
  figmessage("Figure 2G: m.p amp vs. distance")
  vrest <- -70
  plot.feature.dist(r, feature=r$vsrimax_mean - vrest, liw=r$vsrimax_stderr,
                    ylab=expression("Peak m.p. (mV)"))
  panlabel("G", line=-0.3)

  figmessage("Figure 2H: m.p integral vs. distance")
  plot.feature.dist(r, feature=r$vsriint_mean-r$tstop*vrest,
                    liw=r$vsriint_stderr,
                    sem.max=10,
                    ylab=expression(paste("m.p. integral (", mu, "Vs)")))
  panlabel("H", line=-0.3)
  
  figmessage("Figure 2I: m.p. delay vs. distance")
  plot.feature.dist(r, feature=r$vsridel_mean, liw=r$vsridel_stderr,
                    ylab=expression("Delay to peak m.p. (ms)"),
                    legend.pos="topleft")
  panlabel("I", line=-0.3)
  
  ## Bottom row: delay

  figmessage("Figure 2J: Ca amp vs. distance")
  plot.feature.dist(r, feature=r$casrimax_mean*1000, liw=r$casrimax_stderr,
                    ylab=expression(paste("Peak [Ca] (", mu, "M)",sep="")))
  panlabel("J", line=-0.3)

  figmessage("Figure 2K: Ca integral vs. distance")
  plot.feature.dist(r, feature=r$casriint_mean, liw=r$casriint_stderr,
                    ylab=expression(paste("[Ca] integral (", mu, "Ms)",sep="")))
  panlabel("K", line=-0.3)
  
  figmessage("Figure 2L: Ca delay vs. distance")
  plot.feature.dist(r, feature=r$casridel_mean, liw=r$casridel_stderr,
                    ylab=expression("Delay to peak [Ca] (ms)"),
                    legend.pos="topleft")
  panlabel("L", line=-0.3)
  
  dev.off()


  ######################################################################
  ## Figure 4A-H
  ##
  ## Asynchronous stimulation
  ######################################################################
  
  figmessage("Figure 4A-H: Measures of distance with asynchronous SC stimulation")
  
  dataset <- "ca1_poirazi-dendburst-s240-j10t1-a200-n45-bv-r170-sc0-Ra050-nr0100-cv1"
  r <- get.dataset(dataset, dir=file.path("datastore", group))

  stdpostscript(file=dataset, group=group, width=2*colwidth, height=1*colwidth)
  par(mfrow=c(2, 4))
  stdpars()
  par(oma=c(0, 0, 1, 0))
  
  ## 1st row: Asynchronous voltage
  vrest <- -70

  figmessage("Figure 4A: m.p. with asynchronous stimulation")
  plot.synrec.closeup(r$tsoma, r$screc_vsyn,
                      xlab=expression(paste(italic(t), " (ms)")),
                      ylab=expression(paste(italic(V), " (mV)")),
                      panlab="",
                      x.extra=r$tsoma, y.extra=r$vsoma,
                      ylim=c(-80,20),
                      xlim.closeup=NULL)
  panlabel("A", panlab.adj=-0.22)
  title("Asynchronous stimulation", outer=TRUE)

  figmessage("Figure 4B: m.p. amp with asynchronous stimulation")
  plot.feature.dist(r, feature=r$vsrimax_mean - vrest, liw=r$vsrimax_stderr,
                    ylab=expression("Peak m.p. (mV)"))
  panlabel("B")

  figmessage("Figure 4C: m.p. integral with asynchronous stimulation")
  plot.feature.dist(r, feature=r$vsriint_mean-r$tstop*vrest,
                    liw=r$vsriint_stderr,
                    sem.max=10,
                    ylab=expression(paste("m.p. integral (", mu, "Vs)")),
                    legend.pos="topright")
  panlabel("C", panlab.adj=-0.22)

  figmessage("Figure 4D: m.p. amp delay with asynchronous stimulation")
  plot.feature.dist(r, feature=r$vsridel_mean, liw=r$vsridel_stderr,
                    ylab=expression("Delay to peak m.p. (ms)"),
                    legend.pos="topleft")
  panlabel("D")

  ## 2nd row: Asynchronous calcium

  figmessage("Figure 4E: Ca with asynchronous stimulation")
  plot.synrec.closeup(r$tsoma,1000*r$screc_casyn,
                      xlab=expression(paste(italic(t), " (ms)")),
                      ylab=expression(paste("[Ca] (", mu,"M)",sep="")),
                      panlab="",
                      xlim.closeup=NULL,
                      ylim=c(0,100))
  panlabel("E")

  figmessage("Figure 4F: Ca amp with asynchronous stimulation")
  plot.feature.dist(r, feature=r$casrimax_mean*1000, liw=r$casrimax_stderr,
                    ylab=expression(paste("Peak [Ca] (", mu, "M)",sep="")))
  panlabel("F")

  figmessage("Figure 4G: Ca integral with asynchronous stimulation")
  plot.feature.dist(r, feature=r$casriint_mean, liw=r$casriint_stderr,
                    ylab=expression(paste("[Ca] integral (", mu, "Ms)",sep="")),
                    legend.pos="topright")
  panlabel("G")
  
  figmessage("Figure 4H: Ca delay with asynchronous stimulation")
  plot.feature.dist(r, feature=r$casridel_mean, liw=r$casridel_stderr,
                    ylab=expression("Delay to peak [Ca] (ms)"),
                    legend.pos="topleft")
  panlabel("H")
  
  dev.off()

  figmessage("Figure 4I-O: Measures of distance with subthreshold stimulation")
  dataset <-  "ca1_poirazi-dendburst-s190-j00t1-a200-n45-bv-r170-sc0-Ra050-nr0100-cv1"
  r <- get.dataset(dataset, dir=file.path("datastore", group))

  stdpostscript(file=dataset, group=group, width=2*colwidth, height=1*colwidth)
  par(mfrow=c(2, 4))
  stdpars()
  par(oma=c(0, 0, 1, 0))
  
  ## 1st row: Voltage
  vrest <- -70

  figmessage("Figure 4I: m.p. with subthreshold stimulation")
  plot.synrec.closeup(r$tsoma, r$screc_vsyn,
                      xlab=expression(paste(italic(t), " (ms)")),
                      ylab=expression(paste(italic(V), " (mV)")),
                      panlab="",
                      x.extra=r$tsoma, y.extra=r$vsoma,
                      ylim=c(-80,20),
                      xlim.closeup=NULL)
  panlabel("I")
  title("Subthreshold stimulation", outer=TRUE)

  figmessage("Figure 4J: m.p. amp with subthreshold stimulation")
  plot.feature.dist(r, feature=r$vsrimax_mean - vrest, liw=r$vsrimax_stderr,
                    ylab=expression("Peak m.p. (mV)"))
  panlabel("J")

  figmessage("Figure 4K: m.p. amp delay with subthreshold stimulation")
  plot.feature.dist(r, feature=r$vsriint_mean-r$tstop*vrest,
                    liw=r$vsriint_stderr,
                    sem.max=10,
                    ylab=expression(paste("m.p. integral (", mu, "Vs)")),
                    legend.pos="topright")
  panlabel("K")

  figmessage("Figure 4L: m.p. amp delay with subthreshold stimulation")
  plot.feature.dist(r, feature=r$vsridel_mean, liw=r$vsridel_stderr,
                    ylab=expression("Delay to peak m.p. (ms)"),
                    legend.pos="topright")
  panlabel("L")

  
  ## 2nd row: Calcium

  figmessage("Figure 4M: Ca with subthreshold stimulation")
  plot.synrec.closeup(r$tsoma,1000*r$screc_casyn,
                      xlab=expression(paste(italic(t), " (ms)")),
                      ylab=expression(paste("[Ca] (", mu,"M)",sep="")),
                      panlab="",
                      xlim.closeup=NULL,
                      ylim=c(0,100))
  panlabel("M")

  figmessage("Figure 4N: Ca amp with subthreshold stimulation")
  plot.feature.dist(r, feature=r$casrimax_mean*1000, liw=r$casrimax_stderr,
                    ylab=expression(paste("Peak [Ca] (", mu, "M)",sep="")))
  panlabel("N")


  figmessage("Figure 4O: Ca delay with subtrheshold stimulation")
  plot.feature.dist(r, feature=r$casriint_mean, liw=r$casriint_stderr,
                    ylab=expression(paste("[Ca] integral (", mu, "Ms)",sep="")),
                    legend.pos="topright")
  panlabel("O")

  figmessage("Figure 4P: Ca delay with subtrheshold stimulation")
  plot.feature.dist(r, feature=r$casridel_mean, liw=r$casridel_stderr,
                    ylab=expression("Delay to peak [Ca] (ms)"),
                    legend.pos="topright")
  panlabel("P")

  dev.off()
  
  ######################################################################
  ## Figure 6
  ##
  ## Attenuation  measures
  ######################################################################

  figmessage("Figure 6A: Attenuation versus distance")

  dataset <- "epspamp-Ra050"
  r <- get.dataset(dataset, dir=file.path("datastore", group))
  att <- r$epsp_attenuation

  dataset <-  "ca1_poirazi-dendburst-s240-j00t1-a200-n45-bv-r170-sc0-Ra050-nr0100-cv1"
  r <- get.dataset(dataset, dir=file.path("datastore", group))
  r$epsp_attenuation <- att
  
  stdpostscript(file="att-dist", group=group,
                width=0.5*colwidth, height=0.8*colwidth)
  stdpars()

  ## att <- 1-(r$transimp/r$transimp[1])
  
  plot.feature(r, x=att, feature=r$distances,
               ylab=expression(paste("Path distance (", mu, "m)", sep="")),
               ylim=c(0, max(r$distances)), 
               xlab="Attenuation", xlim=c(0, 0.8),
               fit=FALSE,
               plotinds=1:length(att))
  panlabel("A")
  dev.off()

  ## Figure 6 - the same as figure 2, but plotted against attenuation
  figmessage("Figure 6B-E: Synchronous stimulation: Attenuation measures")
  
  stdpostscript(file=paste(dataset, "-att", sep=""), group=group,
                width=1.5*colwidth, height=0.8*colwidth)
  par(mfrow=c(2, 3))
  stdpars()

  ## Top row: peak
  figmessage("Figure 6B: m.p. delay vs. attenuation")
  vrest <- -70
  plot.feature.att(r, feature=r$vsrimax_mean - vrest, liw=r$vsrimax_stderr,
                    ylab=expression("Peak m.p. (mV)"))
  panlabel("B")

  figmessage("Figure 6C: m.p. integral vs. attenuation")
  plot.feature.att(r, feature=r$vsriint_mean-r$tstop*vrest,
                   liw=r$vsriint_stderr,
                   sem.max=10,
                   ylab=expression(paste("m.p. integral (", mu, "Vs)")),
                   legend.pos="topright")
  panlabel("C")
  
  figmessage("Figure 6D: m.p. delay vs. attenuation")
  plot.feature.att(r, feature=r$vsridel_mean, liw=r$vsridel_stderr,
                   ylab=expression("Delay to peak m.p. (ms)"),
                   legend.pos="topleft")
  panlabel("D")

  figmessage("Figure 6E: Ca amplitude vs. attenuation")
  plot.feature.att(r, feature=r$casrimax_mean*1000, liw=r$casrimax_stderr,
                    ylab=expression(paste("Peak [Ca] (", mu, "M)",sep="")))
  panlabel("E")

  figmessage("Figure 6F: Ca integral vs. attenuation")
  plot.feature.att(r, feature=r$casriint_mean, liw=r$casriint_stderr,
                   ylab=expression(paste("[Ca] integral (", mu, "Ms)",sep="")),
                   legend.pos="topright")
  panlabel("F")
  
  figmessage("Figure 6G: Ca delay vs. attenuation")
  plot.feature.att(r, feature=r$casridel_mean, liw=r$casridel_stderr,
                   ylab=expression("Delay to peak [Ca] (ms)"),
                   legend.pos="topleft")
  panlabel("G")
  dev.off()

  ######################################################################
  ##
  ## Figure 7 - Scaled synapses
  ##
  ######################################################################

  dataset <-  "ca1_poirazi-dendburst-s240-j00t1-a200-n45-bv-r170-sc1-Ra050-nr0100-cv1"
  r <- get.dataset(dataset, dir=file.path("datastore", group))
  ## Save for KS test later
  r.sc1 <- r
  
  figmessage("Figure 7A-F: Scaled synapses")
  stdpostscript(file=dataset, group=group, width=colwidth*1.5, height=3)
  par(mfrow=c(2, 3))
  stdpars()

  sem.max.sc <- 5
  ## Top row: peak
  figmessage("Figure 7A: m.p amp vs. distance")
  vrest <- -70
  plot.feature.dist(r, feature=r$vsrimax_mean - vrest, liw=r$vsrimax_stderr,
                    ylab=expression("Peak m.p. (mV)"), sem.max=sem.max.sc)
  panlabel("A", line=-0.4)

  figmessage("Figure 7B: m.p. integral vs. distance")
  plot.feature.dist(r, feature=r$vsriint_mean-r$tstop*vrest,
                    liw=r$vsriint_stderr,
                    sem.max=10*sem.max.sc,
                    ylab=expression(paste("m.p. integral (", mu, "Vs)")),
                    legend.pos="topright")
  panlabel("B", line=-0.4)
  
  figmessage("Figure 7C: m.p. delay vs. distance")
  plot.feature.dist(r, feature=r$vsridel_mean, liw=r$vsridel_stderr,
                    ylab=expression("Delay to peak m.p. (ms)"),
                    legend.pos="topright", sem.max=sem.max.sc)
  panlabel("C", line=-0.4)
  
  ## Bottom row: delay

  figmessage("Figure 7D: Ca amp vs. distance")
  plot.feature.dist(r, feature=r$casrimax_mean*1000, liw=r$casrimax_stderr,
                    ylab=expression(paste("Peak [Ca] (", mu, "M)",sep="")),
                    sem.max=sem.max.sc)
  panlabel("D", line=-0.4)

  figmessage("Figure 7E: Ca integral vs. distance")
  plot.feature.dist(r, feature=r$casriint_mean*1000, liw=r$casriint_stderr,
                    ylab=expression(paste("[Ca] integral (", mu, "Ms)",sep="")),
                    sem.max=sem.max.sc)
  panlabel("E", line=-0.4)
  
  figmessage("Figure 7F: Ca delay vs. distance")
  plot.feature.dist(r, feature=r$casridel_mean, liw=r$casridel_stderr,
                    ylab=expression("Delay to peak [Ca] (ms)"),
                    legend.pos="topright", sem.max=sem.max.sc)
  panlabel("F", line=-0.4)
  
  dev.off()

  ######################################################################
  ##
  ## Figure S1: Comparison of various stimulation protocols
  ##
  ######################################################################

  stdpostscript(file=paste("stimulation-comparison", sep=""), group=group,
                width=1.5*colwidth, height=2*colwidth)
  par(mfcol=c(4, 3))
  stdpars()
  
  figmessage("Figure S1B-E: Dendritic stimulation")
  dataset <-  "ca1_poirazi-dendburst-s240-j00t1-a200-n45-bv-r170-sc0-Ra050-nr0100-cv1"
  r <- get.dataset(dataset, dir=file.path("datastore", group))
  
  figmessage("Figure S1B: m.p amp vs. distance")
  vrest <- -70
  plot.feature.dist(r, feature=r$vsrimax_mean - vrest, liw=r$vsrimax_stderr,
                    ylim=c(0, 100),
                    ylab=expression("Peak m.p. (mV)"))
  panlabel("B")
  
  figmessage("Figure S1C: m.p. delay vs. distance")
  plot.feature.dist(r, feature=r$vsridel_mean, liw=r$vsridel_stderr,
                    ylab=expression("Delay to peak m.p. (ms)"),
                    legend.pos="topleft")
  panlabel("C")

  figmessage("Figure S1D: Ca amp vs. distance")
  plot.feature.dist(r, feature=r$casrimax_mean*1000, liw=r$casrimax_stderr,
                    ylim=c(0,100),
                    ylab=expression(paste("Peak [Ca] (", mu, "M)",sep="")))
  panlabel("D")
  
  figmessage("Figure S1E: Ca delay vs. distance")
  plot.feature.dist(r, feature=r$casridel_mean, liw=r$casridel_stderr,
                    ylab=expression("Delay to peak [Ca] (ms)"),
                    legend.pos="topleft")
  panlabel("E")

  figmessage("Figure S1G-J: Dendritic + somatic stimulation")
  dataset <- "ca1_poirazi-dendburst-s005-j00t1-a200-n45-bv-r170-sc0-Ra050-nr0500-cv1-id63"

  r <- get.dataset(dataset, dir=file.path("datastore", group))

  figmessage("Figure S1G: m.p amp vs. distance")
  vrest <- -70
  plot.feature.dist(r, feature=r$vsrimax_mean - vrest, liw=r$vsrimax_stderr,
                    ylim=c(0, 100),
                    ylab=expression("Peak m.p. (mV)"))
  panlabel("G")

  figmessage("Figure S1H: m.p. delay vs. distance")
  plot.feature.dist(r, feature=r$vsridel_mean, liw=r$vsridel_stderr,
                    ylab=expression("Delay to peak m.p. (ms)"),
                    legend.pos="topleft")
  panlabel("H")

  figmessage("Figure S1I: Ca amp vs. distance")
  plot.feature.dist(r, feature=r$casrimax_mean*1000, liw=r$casrimax_stderr,
                    ylim=c(0,100),
                    ylab=expression(paste("Peak [Ca] (", mu, "M)",sep="")))
  panlabel("I")

  figmessage("Figure S1J: Ca delay vs. distance")
  plot.feature.dist(r, feature=r$casridel_mean, liw=r$casridel_stderr,
                    ylab=expression("Delay to peak [Ca] (ms)"),
                    legend.pos="topleft")
  panlabel("J")
  
  figmessage("Figure S1L-M: Somatic stimulation")

  dataset <- "ca1_poirazi-dendburst-s190-j00t1-a00-n00-bv-r170-sc0-Ra050-nr0001-cv1"
  r <- get.dataset(dataset, dir=file.path("datastore", group))
  
  figmessage("Figure S1L: m.p amp vs. distance")
  vrest <- -70
  plot.feature.dist(r, feature=r$vtreemax_mean-vrest, liw=r$vtreemax_stderr,
                    n.rej=1,
                    ylim=c(0, 100),
                    plotinds=1:length(r$vtreemax_mean), 
                    ylab=expression("Peak m.p. (mV)"))
  panlabel("L")
  
  ## figmessage("Figure S1M: m.p. delay vs. distance")
  plot.new()
  ## plot.feature.dist(r, feature=r$vsridel_mean, liw=r$vsridel_stderr,
  ##                   ylab=expression("Delay to peak m.p. (ms)"),
  ##                   legend.pos="topleft")
  ## panlabel("M")

  figmessage("Figure S1M: Ca amp vs. distance")
  plot.feature.dist(r, feature=r$casrimax_mean*1000, liw=r$casrimax_stderr,
                    n.rej=1,
                    ylim=c(0,100),
                    ylab=expression(paste("Peak [Ca] (", mu, "M)",sep="")))
  panlabel("M")
  
  ## figmessage("Figure S1O: Ca delay vs. distance")
  plot.new()
  ## plot.feature.dist(r, feature=r$casridel_mean, liw=r$casridel_stderr,
  ##                   ylab=expression("Delay to peak [Ca] (ms)"),
  ##                   legend.pos="topleft")
  ## panlabel("O")
  dev.off()  

  ## This figure isn't possible becase the special synapses aren't defined.
  ## figmessage("Figure S1 extra: Synrecs for  dendritic + somatic stimulation")
  ## dataset <- "ca1_poirazi-dendburst-s005-j00t1-a200-n45-bv-r170-sc0-Ra050-nr0100-cv1"
  ## r <- get.dataset(dataset, dir=file.path("datastore", group))
  ## stdpostscript(file=paste(dataset, "-synrecs", sep=""), group=group,
  ##               width=1*colwidth, height=3.5)
  ## plot.synrecs(r, panlabs=c("B", "C", "D", "E", "F"),
  ##              show.IR=TRUE)
  ## dev.off()
  
  ######################################################################
  ##
  ## Figure S3
  ##
  ######################################################################
  
  dataset <-  "ca1_poirazi-dendburst-s240-j00t1-a200-n45-bv-r000-sc0-Ra050-nr0100-cv1"
  r <- get.dataset(dataset, dir=file.path("datastore", group))

  ## Figure S3 - time courses from spines
  figmessage("Figure S3B-F: Time courses in spines")
  stdpostscript(file=paste(dataset, "-synrecs", sep=""), group=group,
                width=colwidth, height=4)
  par(mar=c(1.5, 2.5, 0.1, 0.8))
  plot.synrecs(r, panlabs=c("B", "C", "D", "E", "F"), xlim.closeup=c(28, 34),
               show.IR=TRUE)
  dev.off()

  ## Figure S3 - Measures of distance with SC stimulation
  figmessage("Figure S3G-L: Measures of distance with SC stimulation")
  stdpostscript(file=dataset, group=group, width=1.5*colwidth, height=3)
  par(mfrow=c(2, 3))
  stdpars()
  
  ## Top row: peak
  figmessage("Figure S3G: m.p amp vs. distance")
  vrest <- -70
  plot.feature.dist(r, feature=r$vsrimax_mean - vrest, liw=r$vsrimax_stderr,
                    ylab=expression("Peak m.p. (mV)"))
  panlabel("G", line=-0.3)

  figmessage("Figure S3H: m.p integral vs. distance")
  plot.feature.dist(r, feature=r$vsriint_mean-r$tstop*vrest,
                    liw=r$vsriint_stderr,
                    sem.max=10,
                    ylab=expression(paste("m.p. integral (", mu, "Vs)")))
  panlabel("H", line=-0.3)
  
  figmessage("Figure S3I: m.p. delay vs. distance")
  plot.feature.dist(r, feature=r$vsridel_mean, liw=r$vsridel_stderr,
                    ylab=expression("Delay to peak m.p. (ms)"),
                    legend.pos="topleft")
  panlabel("I", line=-0.3)
  
  ## Bottom row: delay

  figmessage("Figure S3J: Ca amp vs. distance")
  plot.feature.dist(r, feature=r$casrimax_mean*1000, liw=r$casrimax_stderr,
                    ylab=expression(paste("Peak [Ca] (", mu, "M)",sep="")))
  panlabel("J", line=-0.3)

  figmessage("Figure S3K: Ca integral vs. distance")
  plot.feature.dist(r, feature=r$casriint_mean, liw=r$casriint_stderr,
                    ylab=expression(paste("[Ca] integral (", mu, "Ms)",sep="")))
  panlabel("K", line=-0.3)
  
  figmessage("Figure S3L: Ca delay vs. distance")
  plot.feature.dist(r, feature=r$casridel_mean, liw=r$casridel_stderr,
                    ylab=expression("Delay to peak [Ca] (ms)"),
                    legend.pos="topleft")
  panlabel("L", line=-0.3)
  
  dev.off()



  figmessage("Figure 4 extra: Measures of distance with subthreshold stimulation and scaled synapses")
  dataset <-  "ca1_poirazi-dendburst-s190-j00t1-a200-n45-bv-r170-sc1-Ra050-nr0100-cv1"
  r <- get.dataset(dataset, dir=file.path("datastore", group))

  stdpostscript(file=dataset, group=group, width=1.5*colwidth, height=0.75*colwidth)
  par(mfrow=c(2, 4))
  stdpars()
  par(oma=c(0, 0, 1, 0))
  
  ## 1st row: Asynchronous voltage
  vrest <- -70

  figmessage("Figure 4A: m.p. with asynchronous stimulation")
  plot.synrec.closeup(r$tsoma, r$screc_vsyn,
                      xlab=expression(paste(italic(t), " (ms)")),
                      ylab=expression(paste(italic(V), " (mV)")),
                      panlab="",
                      x.extra=r$tsoma, y.extra=r$vsoma,
                      ylim=c(-80,20),
                      xlim.closeup=NULL)
  panlabel("A", panlab.adj=-0.22)
  title("Subtrheshold stimulation with scaled synapses", outer=TRUE)

  figmessage("Figure 4B: m.p. amp with asynchronous stimulation")
  plot.feature.dist(r, feature=r$vsrimax_mean - vrest, liw=r$vsrimax_stderr,
                    ylab=expression("Peak m.p. (mV)"))
  panlabel("B")

  figmessage("Figure 4C: m.p. amp delay with asynchronous stimulation")
  plot.feature.dist(r, feature=r$vsridel_mean, liw=r$vsridel_stderr,
                    ylab=expression("Delay to peak m.p. (ms)"),
                    legend.pos="topleft")
  panlabel("C")

  figmessage("Figure 4D: m.p. integral with asynchronous stimulation")
  plot.feature.dist(r, feature=r$vtreeint_mean-r$tstop*vrest,
                    liw=r$vtreeint_stderr,
                    ylab=expression(paste("m.p. integral (", mu, "Vs)")),
                    legend.pos="topright")
  panlabel("D", panlab.adj=-0.22)

  
  ## 2nd row: Asynchronous calcium

  figmessage("Figure 4E: Ca with asynchronous stimulation")
  plot.synrec.closeup(r$tsoma,1000*r$screc_casyn,
                      xlab=expression(paste(italic(t), " (ms)")),
                      ylab=expression(paste("[Ca] (", mu,"M)",sep="")),
                      panlab="",
                      xlim.closeup=NULL,
                      ylim=c(0,100))
  panlabel("E")

  figmessage("Figure 4F: Ca amp with asynchronous stimulation")
  plot.feature.dist(r, feature=r$casrimax_mean*1000, liw=r$casrimax_stderr,
                    ylab=expression(paste("Peak [Ca] (", mu, "M)",sep="")))
  panlabel("F")

  figmessage("Figure 4G: Ca delay with asynchronous stimulation")
  plot.feature.dist(r, feature=r$casridel_mean, liw=r$casridel_stderr,
                    ylab=expression("Delay to peak [Ca] (ms)"),
                    legend.pos="topleft")
  panlabel("G")

  figmessage("Figure 4H: Ca integral with asynchronous stimulation")
  plot.feature.dist(r, feature=r$casriint_mean, liw=r$casriint_stderr,
                    ylab=expression(paste("[Ca] integral (", mu, "Ms)",sep="")),
                    legend.pos="topright")
  panlabel("H")

  dev.off()


  ######################################################################
  ##
  ## KS test of scaled and unscaled distributions
  ##
  ######################################################################
  figmessage("KS test of scaled and unscaled distributions")
  print(ks.test(r.sc0$casrimax_mean, r.sc1$casrimax_mean))
  
}

## Function to export data to CSV format
export.data <- function(group="car-mag") {
  dataset <-  "ca1_poirazi-dendburst-s240-j00t1-a200-n45-bv-r170-sc0-Ra050-nr0100-cv1"
  r <- get.dataset(dataset, dir=file.path("datastore", group))
  vrest <- -70

  plotinds <- which(r$vsridel_mean!=0)
  fig2h <- cbind(dist=r$dist[plotinds], vsridel_mean=r$vsridel_mean[plotinds])
  write.csv(fig2h, "fig2h.csv", row.names=FALSE)

    plotinds <- which(r$casridel_mean!=0)
  fig2i <- cbind(dist=r$dist[plotinds], casridel_mean=r$casridel_mean[plotinds])
  write.csv(fig2i, "fig2i.csv", row.names=FALSE)
}

plot.all.features.dist.size(group="")
