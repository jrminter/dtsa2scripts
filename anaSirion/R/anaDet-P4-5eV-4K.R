# anaDet-P4-5eV-4K.R
rm(list=ls())

library(qcc)

gitDir <- Sys.getenv("GIT_HOME")
rDir <- paste0(gitDir, '/dtsa2Scripts/anaSirion/R/')
setwd(rDir)

csvFil <- '../dat/Oxford-P4-05eV-4K-15kV-Cu.csv'

df <- read.csv(csvFil,header=TRUE, as.is=TRUE)
print(head(df))

date <- df$date

cw <- df$cw.ev.mu
zo <- df$zo.ev.mu
mn <- df$mn.res.mu
cu.la <- df$cu.la.cts.per.na.sec.mu

qcw <- qcc(cw, type='xbar.one', labels=date,data.name="[eV/ch]")

qzo <- qcc(zo, type='xbar.one', labels=date, data.name="Zero Offset [eV]")

qmn <- qcc(mn, type='xbar.one', labels=date, data.name="Mn Resolution [eV]")

qculA <- qcc(cu.la, type='xbar.one', labels=date, data.name="Cu-La [cts/sec-nA]")


