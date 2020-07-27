# anaDet-P4-5eV-2K.R
rm(list=ls())

library(qcc)

gitDir <- Sys.getenv("GIT_HOME")
rDir <- paste0(gitDir, '/dtsa2Scripts/anaSirion/R/')
setwd(rDir)

csvFil <- '../dat/Oxford-P4-05eV-2K-15kV-Cu.csv'

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

print(date)
df$date <- as.Date(date, "%Y-%m-%d")
print(df)


require(ggplot2)
pd <- position_dodge(0.1) 
dcw <- ggplot( data = df, aes( date, cw.ev.mu )) 
dcw <- dcw + geom_errorbar(aes(ymin=cw.ev.mu-cw.ev.unc,
                               ymax=cw.ev.mu+cw.ev.unc),
                           width=.5, size=1, position=pd)
dcw <- dcw + geom_line(size=1)
dcw <- dcw + geom_point(size=2) 
dcw <- dcw + xlab("Date")
dcw <- dcw + ylab("Channel Width [eV]")
dcw <- dcw + ggtitle("Oxford-P4-05eV-2K")
print(dcw)




dzo <- ggplot(data = df, aes(date, zo.ev.mu)) 
dzo <- dzo + geom_errorbar(aes(ymin=zo.ev.mu-zo.ev.unc,
                               ymax=zo.ev.mu+zo.ev.unc),
                           width=.5, size=1, position=pd)
dzo <- dzo + geom_line(size=1)
dzo <- dzo + geom_point(size=2) 
dzo <- dzo + xlab("Date")
dzo <- dzo + ylab("Zero Offset [eV]")
dzo <- dzo + ggtitle("Oxford-P4-05eV-2K")
print(dzo)


dmn <- ggplot(data = df, aes(date, mn.res.mu)) 
dmn <- dmn + geom_errorbar(aes(ymin=mn.res.mu-mn.res.unc,
                               ymax=mn.res.mu+mn.res.unc),
                           width=.5, size=1, position=pd)
dmn <- dmn + geom_line(size=1)
dmn <- dmn + geom_point(size=2) 
dmn <- dmn + xlab("Date")
dmn <- dmn + ylab("Resolution at MnKa [eV]")
dmn <- dmn + ggtitle("Oxford-P4-05eV-2K")
print(dmn)


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
culacpsna <- culacpsna + ggtitle("Oxford-P4-05eV-2K")
print(culacpsna)



