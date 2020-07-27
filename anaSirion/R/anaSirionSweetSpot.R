# anaSirionSweetSpot.R
rm(list=ls())

gitDir <- Sys.getenv("GIT_HOME")
rDir <- paste0(gitDir, '/dtsa2Scripts/anaSirion/R/')
setwd(rDir)

csvFil <- '../dat/2016-06-17-30kV-Mn-HR-det-posn-analysis.csv'

df <- read.csv(csvFil,header=TRUE, as.is=TRUE)
print(head(df))
plot(mn.ka.cts.per.na.sec.mu~dist.mm, data=df, type='n',
     xlab='working distance [mm]',
     ylab=expression(paste("MnK", alpha, " counts per nA sec")))
points(mn.ka.cts.per.na.sec.mu~dist.mm, data=df, pch=19)
arrows(df$dist.mm,
       df$ mn.ka.cts.per.na.sec.mu - df$mn.ka.cts.per.na.sec.unc, 
       df$dist.mm,
       df$mn.ka.cts.per.na.sec.mu + df$mn.ka.cts.per.na.sec.unc,
       length=0.05, angle=90, code=3)
