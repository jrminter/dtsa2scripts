rm(list=ls())
library(rEDS)

eds.det.panel.plot <- function(csvPath,
                               detTitle,
                               off.cw=0.00015,
                               off.zo=0.20,
                               off.mn=0.15,
                               off.cula=200.,
                               lw.bars=2,
                               len.bars=0.075){
  o.mfrow <- par("mfrow")
  o.mar <- par("mar")
  o.ma <- par("oma")
  par(mfrow=c(2, 2))
  par(mar=c(3.1,4.1,1.0,0.5))
  par(oma=c(0,0,2,0))
  df <- read.csv(csvPath,header=TRUE, as.is=TRUE)
  date <- df$date
  d.min <- min(date)
  d.max <- max(date)
  y.l <- df$cw.ev.mu - df$cw.ev.unc
  y.u <- df$cw.ev.mu + df$cw.ev.unc
  df$date <- as.Date(date, "%Y-%m-%d")
  lims <- data.frame(date=c(min(df$date), max(df$date)),
                     y=c(min(y.l)-off.cw, max(y.u)+off.cw))
  plot(y~date, data=lims, type='n',
       xlab='date', ylab='channel width [ev]')
  points(cw.ev.mu~date, data=df, pch=19)
  arrows(df$date, df$cw.ev.mu,
         df$date, y.l, angle=90, code=2,
         length = len.bars, lwd = lw.bars)
  arrows(df$date, df$cw.ev.mu,
         df$date, y.u, angle=90, code=2,
         length = len.bars, lwd = lw.bars)
  abline(h=mean(df$cw.ev.mu), lw=2, col='blue')
  muCW <- mean(df$cw.ev.mu)
  sMuCW <- sprintf("mean: %.5f", muCW)
  legend("topleft", sMuCW, lw=2, col='blue')

  y.l <- df$zo.ev.mu - df$zo.ev.unc
  y.u <- df$zo.ev.mu + df$zo.ev.unc
  lims <- data.frame(date=c(min(df$date), max(df$date)),
                     y=c(min(y.l)-off.zo, max(y.u)+off.zo))
  plot(y~date, data=lims, type='n',
       xlab='date', ylab='zero offset [ev]')
  points(zo.ev.mu~date, data=df, pch=19)
  arrows(df$date, df$zo.ev.mu,
         df$date, y.l, angle=90, code=2,
         length = len.bars, lwd = lw.bars)
  arrows(df$date, df$zo.ev.mu,
         df$date, y.u, angle=90, code=2,
         length = len.bars, lwd = lw.bars)
  abline(h=mean(df$zo.ev.mu), lw=2, col='blue')
  muZO <- mean(df$zo.ev.mu)
  sMuZO <- sprintf("mean: %.3f", muZO)
  legend("topleft", sMuZO, lw=2, col='blue')

  y.l <- df$mn.res.mu - df$mn.res.unc
  y.u <- df$mn.res.mu + df$mn.res.unc
  lims <- data.frame(date=c(min(df$date), max(df$date)),
                     y=c(min(y.l)-off.mn, max(y.u)+off.mn))
  plot(y~date, data=lims, type='n',
       xlab='date', ylab='FWHM at MnKa [ev]')
  points(mn.res.mu~date, data=df, pch=19)
  arrows(df$date, df$mn.res.mu,
         df$date, y.l, angle=90, code=2,
         length = len.bars, lwd = lw.bars)
  arrows(df$date, df$mn.res.mu,
         df$date, y.u, angle=90, code=2,
         length = len.bars, lwd = lw.bars)
  abline(h=mean(df$mn.res.mu), lw=2, col='blue')
  muMn <- mean(df$mn.res.mu)
  sMuMn <- sprintf("mean: %.3f", muMn)
  legend("topleft", sMuMn, lw=2, col='blue')

  y.l <- df$cu.la.cts.per.na.sec.mu - df$cu.la.cts.per.na.sec.unc
  y.u <- df$cu.la.cts.per.na.sec.mu + df$cu.la.cts.per.na.sec.unc

  lims <- data.frame(date=c(min(df$date), max(df$date)),
                     y=c(min(y.l)-off.cula, max(y.u)+off.cula))
  plot(y~date, data=lims, type='n',
       xlab='date', ylab='CuLa [cts/nA-sec]')
  points(cu.la.cts.per.na.sec.mu~date, data=df, pch=19)
  arrows(df$date, df$cu.la.cts.per.na.sec.mu,
         df$date, y.l, angle=90, code=2,
         length = len.bars, lwd = lw.bars)
  arrows(df$date, df$cu.la.cts.per.na.sec.mu,
         df$date, y.u, angle=90, code=2,
         length = len.bars, lwd = lw.bars)
  abline(h=mean(df$cu.la.cts.per.na.sec.mu), lw=2, col='blue')
  muCuLa <- mean(df$cu.la.cts.per.na.sec.mu)
  sMuCuLa <- sprintf("mean: %d.", round(muCuLa,0))
  legend("topleft", sMuCuLa, lw=2, col='blue')
  mtext(detTitle, side=3, outer = TRUE, cex = 1.25)

  par(mfrow=o.mfrow)
  par(mar=o.mar)
  par(oma=o.ma)
}

gitDir <- Sys.getenv("GIT_HOME")
rDir <- paste0(gitDir, '/dtsa2Scripts/anaSirion/R/')
setwd(rDir)

csvFil <- '../dat/Oxford-P4-05eV-2K-15kV-Cu.csv'


edsDetectorMultiPlot(csvFil, 'FEI Sirion Oxford-P4-05eV-2K')

old <- dev.cur()
dev.next()
pdf(file="FEI-Sirion-Oxford-P4-05eV-2K.pdf", width=9.0,
    height=6.0, pointsize=12, useDingbats=TRUE)
# use the version from library rEDS
edsDetectorMultiPlot(csvFil, 'FEI Sirion Oxford-P4-05eV-2K')
dev.off()
dev.set(old)

