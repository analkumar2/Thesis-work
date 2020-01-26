source("figures-common.R")

## Set colours of lines
palette(c("black","gray","red","blue","green","orange"))

## Change name to dataset in datastore directory
dataset <- "dendburst-s250-j01t1-a200-n45-bv-r170-sc1-Ra050-nr0005-cv1"
dataset <- "ca1_poirazi-dendburst-s240-j00t1-a200-n45-bv-r170-sc0-Ra050-nr0100-cv1"

## Load the data
r <- get.dataset(dataset, dir=".")

## Uncomment next 3 lines if you want to print to file
colwidth <- 3.6                         # Width of column in paper
stdpostscript(file=dataset, group="", width=colwidth, height=3)
## par(mfrow=c(3, 3))                      # Sets 2x2 grid
stdpars()
  
## Top row: peak
## figmessage("Figure A: m.p amp vs. distance")
## vrest <- -70
## plot.feature.dist(r, feature=r$vsrimax_mean - vrest, liw=r$vsrimax_stderr,
##                   ylab="Peak m.p. amp. (mV)")
## panlabel("A")

## figmessage("Figure B: m.p integral vs. distance")
## plot.feature.dist(r, feature=r$vsriint_mean-r$tstop*vrest,
##                   liw=r$vsriint_stderr,
##                   sem.max=10,
##                   ylab=expression(paste("m.p. integral (", mu, "Vs)")))
## panlabel("B")

## figmessage("Figure C: m.p. delay vs. distance")
## plot.feature.dist(r, feature=r$vsridel_mean, liw=r$vsridel_stderr,
##                   ylab="Delay to peak m.p. (ms)",
##                   legend.pos="topleft")
## panlabel("C")
  
## ## Bottom row: delay

## figmessage("Figure D: Ca amp vs. distance")
## plot.feature.dist(r, feature=r$casrimax_mean*1000, liw=r$casrimax_stderr,
##                   ylab=expression(paste("Peak [Ca] (", mu, "M)",sep="")))
## panlabel("D")

## figmessage("Figure 2E: Ca integral vs. distance")
## plot.feature.dist(r, feature=r$casriint_mean, liw=r$casriint_stderr,
##                   ylab=expression(paste("[Ca] integral (", mu, "Ms)",sep="")))
## panlabel("E")

## figmessage("Figure F: Ca delay vs. distance")
## plot.feature.dist(r, feature=r$casridel_mean, liw=r$casridel_stderr,
##                   ylab="Delay to peak [Ca] (ms)",
##                   legend.pos="topleft")
## panlabel("F")

figmessage("Figure G: m.p width vs. distance")
plot.feature.dist(r, feature=r$vsriwidth_mean,
                  liw=r$vsriwidth_stderr/1000,
                  ylab=expression(paste("m.p. half-width (ms)")))
panlabel("G")

## figmessage("Figure H: m.p width vs. m.p")
## plot(r$vsrimax_mean - vrest, r$vsriwidth_mean,
##      xlab="Peak m.p. amp. (mV)",
##      ylab=expression(paste("m.p. half-width (ms)")))
## panlabel("H")


## Uncomment if printing to file
dev.off()
