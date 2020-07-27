# anaConCu7kV.R
gitDir <- Sys.getenv("GIT_HOME")
print(gitDir)
wd <- paste0(gitDir,'/dtsa2Scripts/sim-C-on-Cu/R')
setwd(wd)

bSavePlots <- FALSE

strCsv <- '../csv/sim-C-on-Cu-7kV.csv'

df <- read.csv(strCsv, header = TRUE, as.is = TRUE)
kCKa.loess <- loess(kCKa~tNm,data=df)
print(kCKa.loess)
kCKa.predict <- predict(kCKa.loess, data.frame(tNm=df$tNm))
kCuLa.loess <- loess(kCuLa~tNm,data=df)
kCuLa.predict <- predict(kCuLa.loess, data.frame(tNm=df$tNm))
print(kCuLa.loess)

do.plot <- function(){
  sca <- 0.7
  oldR <- par("mfrow")
  oldM <- par("mar")
  par(mar=c(4.1, 4.1, 0.5, 0.5))
  par(mfrow=c(2,1))
  plot(kCKa~tNm,data=df, type='n', xlab="C thickness [nm]",
       ylab="C Ka K-ratio", cex.axis=sca, cex.lab=sca)
  points(kCKa~tNm,data=df, pch=19, cex=.3)
  lines(df$tNm, kCKa.predict, col='red', lwd=2)
  text(350,0.15,"Monte Carlo simulation of C on Cu at 7 kV",
       cex=0.8)
  legend("topleft", c("MC3 simulation", "loess smooth"),
         cex=c(.7,.7), pch=c(19, NA), lwd=c(NA,2),
         col=c('black', 'red'))

  plot(kCuLa~tNm,data=df, type='n', xlab="C thickness [nm]",
       ylab="Cu La K-ratio", cex.axis=sca, cex.lab=sca)
  points(kCuLa~tNm,data=df, pch=19, cex=.3)
  lines(df$tNm, kCuLa.predict, col='red', lwd=2)
  text(350,0.9,"Monte Carlo simulation of C on Cu at 7 kV",
       cex=0.8)
  legend("bottomleft", c("MC3 simulation", "loess smooth"),
         cex=c(.7,.7), pch=c(19, NA), lwd=c(NA,2),
         col=c('black', 'red'))

  par(mfrow=oldR)
  par(mar=oldM)
}

c.on.cu.7kV <- data.frame(tNmC=df$tNm, kCKa=kCKa.predict,
                          kCuLa=kCuLa.predict)

strRD <- '../RData/c.on.cu.7kV.RData'
save(c.on.cu.7kV, file=strRD)

do.plot()

if(bSavePlots==TRUE){
  library(rEDS)
  rsPlotPng(do.plot(), "../plt/c-on-cu-7kV.png", width = 7.5,
            height = 5, pts = 16, dpi = 100)
  rsPlotPdf(do.plot(), "../plt/c-on-cu-7kV.pdf", width = 9,
            height = 6, pts = 12)
}
