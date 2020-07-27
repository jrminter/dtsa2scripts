# testPredictThickC.R

rm(list=ls())

gitDir <- Sys.getenv("GIT_HOME")
str.wd <- paste0(gitDir, '/dtsa2Scripts/estThickC/R/')
setwd(str.wd)

predict.c.thick.5kV <- function(c.to.si, debug=FALSE){
  t.C.nm <- seq(from=0.5, to=100, by=0.5)
  df2 <- data.frame(t.C.nm=t.C.nm, c.to.si.mu=t.C.nm)
  fi <- '../dat/rda/c-to-si-loess-5kV.rda'
  load(fi)
  if(debug) print(ls())
  df2$c.to.si.mu <- predict(c.to.si.loess.5kV, newdat=df2)
  df2$dif <- sqrt((df2$c.to.si.mu-c.to.si)^2)
  mv = min(df2$dif)
  t <- df2$t.C.nm[which(df2$dif==mv)]
  return(t)
}

predict.c.thick.7kV <- function(c.to.si, debug=FALSE){
  t.C.nm <- seq(from=0.5, to=100, by=0.5)
  df2 <- data.frame(t.C.nm=t.C.nm, c.to.si.mu=t.C.nm)
  fi <- '../dat/rda/c-to-si-loess-7kV.rda'
  load(fi)
  if(debug) print(ls())
  df2$c.to.si.mu <- predict(c.to.si.loess.7kV, newdat=df2)
  df2$dif <- sqrt((df2$c.to.si.mu-c.to.si)^2)
  mv = min(df2$dif)
  t <- df2$t.C.nm[which(df2$dif==mv)]
  return(t)
}

predict.c.thick.10kV <- function(c.to.si, debug=FALSE){
  t.C.nm <- seq(from=0.5, to=100, by=0.5)
  df2 <- data.frame(t.C.nm=t.C.nm, c.to.si.mu=t.C.nm)
  fi <- '../dat/rda/c-to-si-loess-10kV.rda'
  load(fi)
  if(debug) print(ls())
  df2$c.to.si.mu <- predict(c.to.si.loess.10kV, newdat=df2)
  df2$dif <- sqrt((df2$c.to.si.mu-c.to.si)^2)
  mv = min(df2$dif)
  t <- df2$t.C.nm[which(df2$dif==mv)]
  return(t)
}

predict.c.thick.20kV <- function(c.to.si, debug=FALSE){
  t.C.nm <- seq(from=0.5, to=100, by=0.5)
  df2 <- data.frame(t.C.nm=t.C.nm, c.to.si.mu=t.C.nm)
  fi <- '../dat/rda/c-to-si-loess-20kV.rda'
  load(fi)
  if(debug) print(ls())
  df2$c.to.si.mu <- predict(c.to.si.loess.20kV, newdat=df2)
  df2$dif <- sqrt((df2$c.to.si.mu-c.to.si)^2)
  mv = min(df2$dif)
  t <- df2$t.C.nm[which(df2$dif==mv)]
  return(t)
}


a <- predict.c.thick.5kV(0.5)
print(a)
a <- predict.c.thick.7kV(0.2)
print(a)
a <- predict.c.thick.10kV(0.07)
print(a)
a <- predict.c.thick.20kV(0.0125)
print(a)