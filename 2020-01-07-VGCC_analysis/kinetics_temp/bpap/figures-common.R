if (!require("gplots")) {install.packages("gplots")} 
require("gplots")

###
### Global variables
### 

## Whether to cache datasets. Dangerous if data changes between plots.
dataset.cache <- FALSE

## Set labels
dist.lab <- expression(paste("Path distance (",mu,"m)",sep=""))
dist.lim <- c(0,1000)
size.lab <- "Effective amplitutde (mV)"
size.lim <- c(0,0.25)
att.lab  <- "Attenuation"
att.lim  <- c(0, 0.6)
figdir <- "figures/"

## Standard graphics parameters
stdparlist <-
  list(mar=c(2.2, 2.5, 0.6, 0.8),       # bottom, left, top, right
       mgp=c(1.2, 0.2, 0), # Distance of the axis title, axis labels
                           # and axis line from the edge of the plot
       oma=c(0,0.5,0,0),
       cex=1,                         # This is the match the change
                                      # caused by using mfrow in the
                                      # table plots. It may need to
                                      # change
       tcl=-0.2)                      # Tick length

## Function to set pars of multiple plots
stdpars <- function() {
  par(stdparlist)
  par(cex=1)                         # This is the match the change caused by us
}

## Place a label on a panel
panlabel <- function(panlabel, panlab.adj=-0.18, line=-1) {
  ##
  ## mtext(text=panlabel, font=2, line=-1, side=3, adj=panlab.adj, cex=8/6)
  mtext(panlabel, side=3, adj=-par("plt")[1]/(par("plt")[2]-par("plt")[1]),
        line=line, font=2, cex=8/8)

  ## mtext(text=panlabel,font=2,line=2.75,side=2,padj=panlabel.adj,las=2)
  ## mtext(text=panlabel, font=2,
  ## line=2.75, side=2, adj=1, padj=panlabel.adj, las=0) 
}

## Get a dataset in "R" format, caching it for speed
get.dataset <- function(dataset.str, r=NULL, dir="datastore") {
  dataset <- gsub("-",".",dataset.str) # need to get rid of - for R
  if (!exists(dataset) || !dataset.cache) {
    file <- file.path(dir, paste(dataset.str,".R", sep=""))
    print(paste("Reading",file))
    source(file ,local=TRUE)
    eval(parse(text=paste(dataset,"<<-r")))
  }
  return(eval(parse(text=dataset)))
}

## Find element of table closest in value to x 
near <- function(x,table) {
  which.min(abs(table-x))
}



epsp.plot <- function(epsp=epsps, colour=TRUE) {
  epsp$vsoma[epsp$vsoma==0] = Erest
  epsp$somaamp <- -Erest + apply(epsp$vsoma,1,max)
  tree.plot(epsp$distance,epsp$somaamp,epsp$diams,epsp$parents,
            xlab = expression(paste("path distance /", mu, "m")) ,
            ylab = expression("EPSP amplitude at soma /mV"), colour=colour)
}

tree.plot.dist.epspsize <- function(epsps=epsps, bpaps,
                                    y=bpaps$synamp, ylab="blee", name,
                                    colour=TRUE,tree=T)
{
  tree.plot(bpaps$distance,y,bpaps$diams,bpaps$parents,
            xlab = expression(paste("path distance /", mu, "m")) ,
            ylab = ylab, colour=colour,tree=tree)
  if (colour==TRUE) {
    filename = paste(name, "-distance-", dataset, ".eps", sep="")
  } else {
    filename = paste(name, "-distance-", dataset, "-bw.eps", sep="")
  }
  dev.copy2eps(file=paste(figdir,filename,sep=""), width=4, height=3)

  tree.plot(epsps$somaamp,y,bpaps$diams,bpaps$parents,
            xlab = "EPSP amplitude at soma /mV",
            ylab = ylab, colour=colour)

  if (colour==TRUE) {
    filename = paste(name, "-epspsomaamp-", dataset, ".eps", sep="")
  } else {
    filename = paste(name, "-epspsomaamp-", dataset, "-bw.eps", sep="")
  }

  dev.copy2eps(file=paste(figdir,filename,sep=""), width=4, height=3)
}

## Get all ancestors of a segments
get.ancestors <- function(r, i) {
  parents <- r$parents + 1
  inds <- c()

  ## Get parents
  while (i != 1) {
    inds <- c(inds, i)
    i <- parents[i]
  }
  inds <- c(inds, i)
  return(inds)
}



## Return a list of the indicies of the parents of the segments
## indexed in the vector SEGS
get.parent.segs <- function(r, segs) {
  parents <- rep(NA, length(segs))
  for (i in 1:length(segs)) {
    parents[i] <- which(get.parent.seg(r,segs[i],segs)==segs)[1]
  }
  return(parents)
}

## Find the parent relationships of segment seg in a vector of segment
## indicies segs
get.parent.seg <- function(r, seg, segs) {
  ## Get the index of the parent segment (+1 converts from NEURON to R indicies)
  parent.seg <- r$parents[seg+1]
  ## Is this segment in the vector of segment indicies?
  if (parent.seg %in% segs) {
    ## If it is, return the index of the parent segment
    return(parent.seg)
  } else {
    ## If not, search. 
    if (parent.seg==0)  {
      ## If the parent segment is 0, this means we are in the root compartment, so don't return anything
      return(NA)
    } else {
      ## Otherwise, see if we can find a parent of the parent that is in the list of segments
      get.parent.seg(r, parent.seg, segs)
    }
  }
}

## Plot a FEATURE (e.g. Maxium calcium concentration or maximum
## membrane potential) versus a quantity X (e.g. distance or EPSP
## amplitude at soma). Upper and lower limits on error bars (UIW and
## LIW) can be set. The dataset R needs to be provided and optionally,
## the curve will be FIT by linear regression. Only datapoints indexed
## by PLOTINDS are plotted. By default, PLOTINDS is set to correspond
## to nonzero points
plot.feature <- function(r, x, feature, liw=NULL, uiw=liw,
                         xlab, ylab, xlim=c(0,1000), ylim=NULL,
                         fit=TRUE, legend.pos="topright",
                         legend.swap=FALSE,
                         plotinds=ifelse(is.null(liw),
                           which(feature != 0),
                           which((feature != 0) & (liw != 0))),
                         verbose=TRUE, epsilon=0.0001,
                         ordinate.dependent=TRUE,
                         sem.max=1,
                         n.rej=10,
                         ...) {
  ## highlights means the "special", coloured synapses
  highlights <- r$screc_inputs_sri+1

  branchcol <- gray(0.4)
  trunkcol <- "black"
  highlightcol <- 2+seq(1:length(highlights))
  cols <- rep(branchcol, length(feature))
  cols[r$trunk_sri+1] <- trunkcol
  cols[highlights] <- highlightcol

  ## Ignore points that are activated less then 10 times (to ensure
  ## reliable mean and SEM)
  if (!is.null(liw)) {
    lpi <- length(plotinds)
    vectorbreaks <- -0.5:(max(r$syn_sri) + 0.5)
    histresults <- hist(r$syn_sri, breaks=vectorbreaks, plot=FALSE)
    plotinds <- plotinds[histresults$counts[plotinds] >= n.rej]
    cat(paste(lpi - length(plotinds), "points rejected because they are activated less than ", n.rej, " times\n")) 
    ## plotinds <- plotinds[liw[plotinds] < sem.max]
    ## cat(paste(lpi - length(plotinds), "points rejected because their SEM is greater than", sem.max, "\n"))
  }
  
  ## Ignore points that are activated less then 10 times (to ensure
  ## reliable mean and SEM)
  
  ## Only plot features selected by plotinds
  feature <- feature[plotinds]
  x <- x[plotinds]
  diams <- r$diams[plotinds]
  parents <- get.parent.segs(r, plotinds-1)
  cols <- cols[plotinds]
  if (!is.null(liw)) {
    liw <- liw[plotinds]
  }

  ## Find correct position for legend
  if (legend.swap) {
    legend.pos.swapped <- switch(legend.pos,
      topleft = "topright",
      topright = "topleft")
  } else {
    legend.pos.swapped <- legend.pos
  }

  if (ordinate.dependent) {
    X <- feature
    Y <- x
    Xlab <- ylab
    Ylab <- xlab
    Xlim <- ylim
    Ylim <- xlim
  } else {
    X <- x
    Y <- feature
    Xlab <- xlab
    Ylab <- ylab
    Xlim <- xlim
    Ylim <- ylim
  }

  ## Do the plot itself
  tree.plot(X, Y, diams, parents,
            xlab=Xlab, ylab=Ylab,
            xlim=Xlim, ylim=Ylim,
            colour=FALSE, col="gray",
            liw=liw, uiw=uiw, err='x',
            highlights=highlights, 
            ...)

  ## Plot the varying classes of points. Do it in a particular order so
  ## that the trunk, branch and highlighted points are distinct.
  i <- which(cols==branchcol)
  points(X[i], Y[i],  col=cols[i], pch=20, cex=0.8)
  i <- which(cols==trunkcol)
  points(X[i], Y[i],  col=cols[i], pch=20, cex=0.8)
  i <- which((cols!=trunkcol) & (cols!=branchcol))
  points(X[i], Y[i], col=cols[i], pch=20, cex=0.8)

  ## Do the curve fitting
  if (fit!=FALSE) {
    ## Compute linear fit
    fmlin <- lm(x ~ y + 1 , data.frame(y=feature, x=x))
    sfmlin <- summary(fmlin)

    ## Compute expontential fit, by starting with log fit...
    fmlog <- lm(log(x+epsilon) ~ y, data.frame(y=feature, x=x))
    sfmlog <- summary(fmlog)
    ## ... and refining with optim
    opt <- optim(fmlog$coef, function(p) { sum(exp(p[1] + p[2]*feature) - x)^2})
    f <- exp(opt$par[1] + opt$par[2]*feature)
    r.squared.log <- 1 - sum((f - x)^2)/sum((x-mean(x))^2)

    cat(paste("Linear R^2 =", format(sfmlin$r.squared, digits=2),
                "; Exponential R^2 = ", format(r.squared.log, digits=2), "\n"))
    
    if (sfmlin$r.squared > r.squared.log) {
      cat("Linear fit selected\n")
      if (ordinate.dependent) {
        abline(fmlin$coef[1], fmlin$coef[2])
      } else {
        abline(-fmlin$coef[1]/fmlin$coef[2], 1/fmlin$coef[2])
      }
      sfm <- sfmlin
      r.squared <- sfm$r.squared
    } else {
      cat("Exponential fit selected\n")
      if (ordinate.dependent) {
        lines(sort(feature), exp(opt$par[1] + opt$par[2]*sort(feature)))
      } else {
        lines(exp(opt$par[1] + opt$par[2]*sort(feature)), sort(feature))
      }
      sfm <- sfmlog
      r.squared <- r.squared.log
    }
    fstatistic <- sfm$fstatistic
    p.value <- pf(fstatistic[1], fstatistic[2], fstatistic[3], lower.tail = FALSE)
    digits <- max(3, getOption("digits") - 3)
    cat("F-statistic:", formatC(fstatistic[1L], digits=digits),
        "on", fstatistic[2L], "and",
        fstatistic[3L], "DF,  p-value:",
        format.pval(p.value, digits=digits),
        "\n")

    if (p.value > 0.001) {
      warning("F statistic not significant at 0.001 level")
    }
    legend(x=legend.pos.swapped,
           legend=c(eval(substitute(expression(paste(italic(R)^2,"=",r2)),list(r2=format(r.squared, digits=2))))),
           adj=c(0.3,0),yjust=0.5)
  }
}

## Plot a FEATURE (e.g. Maxium calcium concentration or maximum
## membrane potential) versus a distance. See plot.feature for
## explanation.
plot.feature.dist <- function(r, feature=r$casrimax_mean*1000,
                              liw=r$casrimax_sterr,
                              plotinds=which((feature != 0) & (liw != 0)),
                              ylab=expression(paste("Peak [Ca] /", mu, "M",sep="")),
                              xlim=NULL,
                              legend.pos="topright", ...) {
  plot.feature(r, x=r$distances, feature=feature, liw=liw,
               ylab=ylab,
               xlim=dist.lim, xlab=dist.lab,
               plotinds=plotinds,
               legend.pos=legend.pos, ...)
}

## Plot a FEATURE (e.g. Maxium calcium concentration or maximum
## membrane potential) versus attenuation. Upper and
## lower limits on error bars (UIW and LIW) can be set. The dataset R
## needs to be provided and optionally, the curve will be FIT by
## linear regression. Any zero datapoints from FEATURE can be removed
## (REMOVE.ZEROS)
plot.feature.att <- function(r, feature=r$casrimax_mean*1000,
                             liw=r$casrimax_sterr,
                             plotinds=which((feature != 0) & (liw != 0)),
                             ylab=expression(paste("Peak [Ca] /", mu, "M",sep="")),
                             xlim=NULL,
                             legend.pos="topright", ...) {
  att <- r$epsp_attenuation
  plot.feature(r, x=att, feature=feature, liw=liw,
               ylab=ylab,
               xlim=att.lim, xlab=att.lab,
               plotinds=plotinds,
               legend.pos=legend.pos, ...)
}


## Core function to plot points with connections between parents and children
tree.plot <- function (x, y, d, parents, liw=NULL, uiw=liw, diam.thresh=1.5,
                       xlab, ylab, colour=TRUE, newaxes=TRUE,
                       col="black", tree=TRUE,  xlim=c(0,1200), err='y',
                       highlights=NULL,ylim=NULL,...) {
  cols <- c()
  if (colour==TRUE) {
    cols[d>diam.thresh]  <- "red"
    cols[d<=diam.thresh] <- "blue"
  } else {
    cols[d>=0] <- "gray"
  }
  ## Point colours
  pcols <- cols
  pcols[highlights] <- seq(1:length(highlights))
  if (newaxes==TRUE) {
    if (is.null(liw)) {
      plot(x, y, col=cols, xlab=xlab, ylab=ylab, pch=".", xlim=xlim, ylim=ylim, ...)
    } else {
      ##      plot(x,y,col=6,xlab=xlab,ylab=ylab,pch=".",xlim=xlim)
      plotCI(x, y, col=cols, xlab=xlab, ylab=ylab,  pch=".",  xlim=xlim, liw=uiw, uiw=uiw, sfrac=0, gap=0, err='x')
    }
  }
  if (tree==TRUE) {
    for (i in 1:length(x)) {
      lines(c(x[i], x[parents[i]]),c(y[i],y[parents[i]]),col=cols[i],lwd=d[i],...)
    }
  }
}

plot.vtree <- function(r,inds) {
  matplot(r$tsoma-40,t(r$vtree[inds,]),type="l",lty=1,lwd=2,
          xlab=expression(paste(italic(t), " (ms)")),ylab="V /mV",xlim=c(0,60),bty="n")
}

plot.casyn <- function(r,inds) {
  matplot(r$tsoma-40,1000*t(r$casyn[inds,]),type="l",lty=1,lwd=2, xlab="t
  /ms",ylab=expression(paste("[Ca] /", mu,
  "M",sep="")),xlim=c(0,60),bty="n")
}

plot.synrec <- function(x, y, xlab=expression(paste(italic(t), " (ms)")), ylab="",
                        tlim=c(40,100),
                        xlim=c(0,tlim[2]-tlim[1]),
                        adj=-0.5,
                        mar=c(4,4,1,1),col=seq(3,dim(y)[2]),lwd=2,
                        x.extra=NULL, y.extra=NULL,lty=1,
                        panlab="", panlab.adj=-0.5, outer=FALSE,
                        ...) {
  oldmar <- par()$mar
  par(mar=mar)
  ## Transform x
  x <- x - tlim[1]
  if (outer) {
    ## Chop out unwanted x and y
    y <- y[,(x >= xlim[1]) & (x <= xlim[2])]
    x <- x[(x >= xlim[1]) & (x <= xlim[2])]
  }
  matplot(x, t(y),type="l",lty=lty,
          lwd=lwd,
          xlab=ifelse(outer, "", xlab), ylab=ylab, xlim=xlim,
          bty="n",col=col, ...)

  if (!is.null(x.extra)) {
    ## Transform x.extra
    x.extra <- x.extra - tlim[1]
    if (outer) {
      ## Chop out unwanted x and y
      y.extra <- y.extra[(x.extra >= xlim[1]) & (x.extra <= xlim[2])]
      x.extra <- x.extra[(x.extra >= xlim[1]) & (x.extra <= xlim[2])]
    }
    lines(x.extra,y.extra,lty=lty,
          lwd=1,
          xlim=xlim,
         col="black")
  }

  if (!is.null(panlab)) {
    panlabel(panlab, panlab.adj)
  }
  ## This does not work
  if (length(xlab) & outer) {
     par(xpd=NA)
     mtext(text=xlab, side=1, line=par("mgp")[1])
     par(xpd=TRUE)
   }
  par(mar=oldmar)
}

## Plot full version of synrec and its close-up
plot.synrec.closeup <- function(x,y,
                                xlab="",xaxt="s",
                                ylab="",
                                panlab="",mar.bottom=par()$mar[1],
                                xlim.closeup=c(28,35), yat=NULL,
                                ...) {
  plot.synrec(x,y,
              xlab=xlab, xaxt=xaxt,
              ylab=ylab,
              panlab=panlab,
              panlab.adj=ifelse(is.null(xlim.closeup),-0.25,-0.25),
              mar=c(mar.bottom, par("mar")[2:4]),
              yaxt=ifelse(!is.null(yat), "n", "s"), ...) #3.5,0.5,1),
  ## if (outer) {
  ##   mtext(xlab, 1, line=0, outer=TRUE)
  ## }
  ## Pretty label
  if (!is.null(yat)) {
    axis(2, at=yat, label=format(yat, zero.print=TRUE))
  }
  if(!is.null(xlim.closeup)) {
    plot.synrec(x,y,
                xlab=xlab, xaxt=xaxt,
                ylab="",yaxt="n",
                mar=c(mar.bottom,0.5, par("mar")[3],1),
                xlim=xlim.closeup,...)
    ## if (outer) {
    ##   mtext(xlab, 1, line=0, outer=TRUE)
    ## }
  }

}

## Plot V, [Ca], I_CaNMDA and I_CaAMPA and (if show.IR is TRUE), I_CaR
## If xlim.closeup is non-null, plot an enlargement of the trace.
plot.synrecs <- function(r, xlim.closeup=c(25,35), panlabs, show.IR=FALSE) {
  ## Number of rows to plot
  nr <- ifelse(show.IR, 5, 4)
  ## Set layout if there is going to be a closeup
  if (!is.null(xlim.closeup)) {
    layout(matrix(1:(nr*2),nr,2,byrow=TRUE),
           widths= c(1,1),
           heights=c(1,1,1,1,1))
  } else {
    par(mfrow=c(nr,1))
  }

  ## Set standard plotting parameters
  stdpars()
  par(mar=c(0.5, par("mar")[2]+1, 0.1, par("mar")[4]))
  par(oma=c(2,0,0,0))

  ## Voltage, V
  plot.synrec.closeup(r$tsoma, r$screc_vsyn,
                      xlab="", xaxt="n",
                      ylab=expression(paste(italic(V), " (mV)")),
                      panlab=panlabs[1],
                      x.extra=r$tsoma, y.extra=r$vsoma,
                      ylim=c(-80,20),
                      xlim.closeup=xlim.closeup, yat=c(-80, -40, 0))

  ## Calcium concentration [Ca]
  plot.synrec.closeup(r$tsoma,1000*r$screc_casyn,
                      xlab="",xaxt="n",
                      ylab=expression(paste("[Ca] (", mu,"M)",sep="")),
                      panlab=panlabs[2],
                      xlim.closeup=xlim.closeup)

  ## I_CaNMDA
  
  ## icasyn is the current density in mA/cm2.  To convert this to
  ## current in nA, it is multiplied by the areas of the spine heads
  ## and an appropriate factor (1E-2 = 1E-8 * 1E6). At the moment this
  ## is not actually plotted; it is here for reference.
  
  area <- r$headdiam * r$headL * pi     # Area of spine head in um^2
  icasyn <- r$screc_icasyn * area * 1E-2

  ica_nmdasyn <- r$screc_ica_nmdasyn
  ica_ampasyn <- r$screc_ica_ampasyn
  ## The factor of 1000 converts nA to pA
  plot.synrec.closeup(r$tsoma,1000*ica_nmdasyn,
                      xlab="",xaxt="n",
                      ylim=c(min(1000*ica_nmdasyn),0),
                      ylab=expression(paste(italic(I)[CaNMDA], " (pA)",sep="")),
                      panlab=panlabs[3],
                      col=c(3,4,5,6),
                      lwd=2,lty=c(1,1,1),
                      xlim.closeup=xlim.closeup, yat=c(0, -0.05))

  ## I_CaAMPA
  plot.synrec.closeup(r$tsoma,1000*r$screc_ica_ampasyn,
                      xlab=ifelse(show.IR, "", expression(paste(italic(t), " (ms)"))),
                      xaxt=ifelse(show.IR, "n", "s"),
                      ylim=c(min(1000*ica_ampasyn),0),
                      ylab=expression(paste(italic(I)[CaAMPA], " (pA)",sep="")),
                      panlab=panlabs[4],
                      col=c(3,4,5,6),
                      lwd=2,lty=c(1,1,1),
                      xlim.closeup=xlim.closeup, yat=c(0, -0.02))

  ## I_CaR
  if (show.IR) {
    icar <- icasyn - ica_ampasyn - ica_nmdasyn
    plot.synrec.closeup(r$tsoma,1000*icar,
                        xlab=expression(paste(italic(t), " (ms)")),
                        ylim=c(min(1000*icar, -0.5),0),
                        ylab=expression(paste(italic(I)[CaR], " (pA)",sep="")),
                        panlab=panlabs[5],
                        col=c(3,4,5,6),
                        lwd=2,lty=c(1,1,1),
                        xlim.closeup=xlim.closeup, yat=c(0, -0.5), outer=TRUE)
  }
  
  return(r$screc_inputs_sri+1)
}

quantities.at.syn.to.quantity.at.seg <- function(r, quant) {
  out <- matrix(NA, nrow(quant), length(r$distances))
  for (i in 1:nrow(quant)) {
    out[i,r$syn_sri[i,] + 1] <- quant[i,]
  }
  return(out)
}

vsyndel.exp <- function() {
  q = quantities.at.syn.to.quantity.at.seg(r, r$vsyndel)
  par(mfcol=c(4, 1))
  stdpars()
  matplot(r$distances, t(q), col='black', pch='.', cex=2, xlab="distance", ylab="Delay to peak voltage")

  mq <- colMeans(q, na.rm=TRUE)
  sdq <- sd(q, na.rm=TRUE)
  Nq <- (apply(!is.na(q), 2, sum))
  print(length(Nq))
  semq <- sdq/(Nq-1)
  print(Nq)
  plotCI(r$distances, mq, uiw=semq, sfrac=0, gap=0, pch=20, ylab="Mean and SEM",
         barcol="blue", xlab="distance")
  points(r$distances[Nq==0], mq[Nq==0], col="red", pch=20)
  legend("topright", c("Mean", "Errorbars show SEM"), text.col=c("black", "blue"))
  plot(r$distances, sd(q, na.rm=TRUE), pch=20, ylab="Standard Deviation", ylim=c(0, 6),
       xlab="distance")
  legend("topright", c("SD", "SEM"), text.col=c("black", "blue"))
  points(r$distances, semq, pch=20, col="blue")
  plot(r$distances, Nq, pch=20, ylab="No. recordings from location",
       xlab="distance")

}

figmessage <- function(mess) {
  s <- "**********************************************************************"
  cat(paste("\n", substr(s, 1, 3), " ",
            mess, " ",
            substr(s, 1, 70 - nchar(mess) - 5), "\n\n",
            sep=""))
}

## Function to set up plot
stdpostscript <- function(file, width=4, height=4, pars=c(), group="car") {
  dir <- file.path("plots", group)
  dir.create(dir, FALSE)
  arialdir <- "/disk/scratch/sterratt/datastore/src/arial"
  Arial <- Type1Font("Arial", c(file.path(arialdir, "arial.afm"),
        file.path(arialdir, "arialbd.afm"),
        file.path(arialdir, "ariali.afm"),
        file.path(arialdir, "arialbi.afm")))
  family <- Arial
  pdf(family=family, file=file.path(dir, paste(file, ".pdf", sep="")),
      width=width, height=height, onefile=FALSE, pointsize=8)

  stdpars()
}
