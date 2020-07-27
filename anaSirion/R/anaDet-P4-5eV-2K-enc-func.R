# anaDet-P4-5eV-2K.R
rm(list=ls())

library(ggplot2)
library(grid)

gitDir <- Sys.getenv("GIT_HOME")
rDir <- paste0(gitDir, '/dtsa2Scripts/anaSirion/R/')
setwd(rDir)

anaDetectorFile <- function(csvPath, detName, bVerbose=FALSE,
                            labOffCW=1, labOffZO=1, labOffMn=1,
                            labOffCu=1){
  df <- read.csv(csvPath,header=TRUE, as.is=TRUE)
  nr <- nrow(df)
  md <- round(nr/2, 0)
  # print(md)
  date <- df$date
  df$date <- as.Date(date, "%Y-%m-%d")
  pd <- position_dodge(0.1)
  v1 <- df$cw.ev.mu
  v2 <- df$cw.ev.unc
  v <- v1+v2
  vMax <- max(v)+ 0.0001
  if(bVerbose){
    print(head(df))
    print(v)
    print(vMax)
  }
  
  dcw <- ggplot( data = df, aes( date, cw.ev.mu )) 
  dcw <- dcw + geom_errorbar(aes(ymin=cw.ev.mu-cw.ev.unc,
                                 ymax=cw.ev.mu+cw.ev.unc),
                             width=.5, size=1, position=pd)
  dcw <- dcw + geom_line(size=1)
  dcw <- dcw + geom_point(size=2) 
  dcw <- dcw + xlab("Date")
  dcw <- dcw + ylab("Channel Width [eV]")
  # dcw <- dcw + ggtitle(detName)
  muCW <- mean(df$cw.ev.mu)
  sMuCW <- sprintf("mean: %.5f", muCW)
  dcw <- dcw + geom_hline(aes(yintercept=muCW),
                          colour='blue', size=1)
  dcw <- dcw + annotate("text", x=df$date[md], y=vMax, # muCW+0.0005,
  # dcw <- dcw + annotate("text", x=df$date[0], y=vMax,
                        label=sMuCW, colour='blue')
  # print(dcw)
  
  v1 <- df$zo.ev.mu
  v2 <- df$zo.ev.unc
  v <- v1+v2
  vMax <- max(v) + .1
  if(bVerbose){
    print(v)
    print(vMax)
  }
  
  dzo <- ggplot(data = df, aes(date, zo.ev.mu)) 
  dzo <- dzo + geom_errorbar(aes(ymin=zo.ev.mu-zo.ev.unc,
                                 ymax=zo.ev.mu+zo.ev.unc),
                             width=.5, size=1, position=pd)
  dzo <- dzo + geom_line(size=1)
  dzo <- dzo + geom_point(size=2) 
  dzo <- dzo + xlab("Date")
  dzo <- dzo + ylab("Zero Offset [eV]")
  # dzo <- dzo + ggtitle(detName)
  muZO <- mean(df$zo.ev.mu)
  sMuZO <- sprintf("mean: %.3f", muZO)
  dzo <- dzo + geom_hline(aes(yintercept=muZO),
                          colour='blue', size=1)
  dzo <- dzo + annotate("text", x=df$date[md], y=vMax, 
                        label=sMuZO, colour='blue')
  # print(dzo)
  
  v1 <- df$mn.res.mu
  v2 <- df$mn.res.unc
  v <- v1+v2
  vMax <- max(v) + .05
  if(bVerbose){
    print(v)
    print(vMax)
  }
  
  dmn <- ggplot(data = df, aes(date, mn.res.mu)) 
  dmn <- dmn + geom_errorbar(aes(ymin=mn.res.mu-mn.res.unc,
                                 ymax=mn.res.mu+mn.res.unc),
                             width=.5, size=1, position=pd)
  dmn <- dmn + geom_line(size=1)
  dmn <- dmn + geom_point(size=2) 
  dmn <- dmn + xlab("Date")
  dmn <- dmn + ylab("Resolution at MnKa [eV]")
  # dmn <- dmn + ggtitle(detName)
  muMn <- mean(df$mn.res.mu)
  sMuMn <- sprintf("mean: %.2f", muMn)
  dmn <- dmn + geom_hline(aes(yintercept=muMn),
                                      colour='blue', size=1)
  dmn <- dmn + annotate("text", x=df$date[md], y=vMax,
                        label =sMuMn, colour='blue')
  # print(dmn)
  
  v1 <- df$cu.la.cts.per.na.sec.mu
  v2 <- df$cu.la.cts.per.na.sec.unc
  v <- v1+v2
  vMax <- max(v) + 100.
  if(bVerbose){
    print(v)
    print(vMax)
  }
 
  culacpsna <- ggplot(data = df, aes(date, cu.la.cts.per.na.sec.mu)) 
  culacpsna <- culacpsna + geom_errorbar(aes(ymin=cu.la.cts.per.na.sec.mu
                                             - cu.la.cts.per.na.sec.unc,
                                             ymax=cu.la.cts.per.na.sec.mu
                                             +cu.la.cts.per.na.sec.unc),
                                         width=.5, size=1, position=pd)
  culacpsna <- culacpsna + geom_line(size=1)
  culacpsna <- culacpsna + geom_point(size=2) 
  culacpsna <- culacpsna + xlab("Date")
  culacpsna <- culacpsna + ylab("Cu-La cps/nA/s")
  # culacpsna <- culacpsna + ggtitle(detName)
  muCuLa <- mean(df$cu.la.cts.per.na.sec.mu)
  sMuCuLa <- sprintf("mean: %.2f", muCuLa)
  culacpsna <- culacpsna + geom_hline(aes(yintercept=muCuLa),
                                          colour='blue', size=1)
  culacpsna <- culacpsna + annotate("text", x=df$date[md], y=vMax,
                                    label=sMuCuLa, colour='blue')
  # print(culacpsna)
  pushViewport(viewport(layout = grid.layout(2, 2)))
  print(dcw, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(dzo, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(dmn, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(culacpsna, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
  
  ret <- list(dcw,dzo, dmn, culacpsna)
  
  return (ret)
}

csvFil <- '../dat/Oxford-P4-05eV-2K-15kV-Cu.csv'
res <- anaDetectorFile(csvFil, "Oxford-P4-05eV-2K",
                       labOffCW=3, labOffZO=3, labOffMn=3, labOffCu=3)

      
      
# print(res) # ['detCW'])

