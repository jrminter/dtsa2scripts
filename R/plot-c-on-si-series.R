# plot-c-on-si-series.R

rm(list=ls())
library(ggplot2)
library(grid)
library(gridExtra)

gitDir <- Sys.getenv("GIT_HOME")
rDir <- paste0(gitDir, '/dtsa2Scripts/R/')
setwd(rDir)


fi <-'../mc3Scripts/dat/csv/C-on-Si-5kV-10000-Traj.csv'

df <- read.csv(fi, header=TRUE, as.is=TRUE)
print(head(df))




df$c.to.si.mu <- df$C.Int.mu / df$Si.Int.mu
df$c.to.si.unc <- sqrt((df$C.Int.unc/df$C.Int.mu)**2 + 
                        (df$Si.Int.unc/df$Si.Int.mu)**2)
df$c.to.si.unc[1]=df$c.to.si.unc[2]

# make a prediction df
t.C.nm <- seq(from=0.5, to=100, by=0.5)
df2 <- data.frame(t.C.nm=t.C.nm)
df3 <- data.frame(t.C.nm=t.C.nm)

siZero <- df$Si.Int.mu[1]
df$Si.Int.rel <- df$Si.Int.mu / siZero

rel.loess <- loess(Si.Int.rel ~ t.C.nm, span=0.75, data=df)
# let's save the fit...
fi <-'../mc3Scripts/dat/rel-loess-5kV.RData'
save(rel.loess, file=fi)

df2$Si.Int.rel = round(predict(rel.loess, newdata = df2), 5)

# let's fix the intercept
intercept <- 1.0
rel.lm <- lm(I(Si.Int.rel - intercept) ~ 0 + t.C.nm, data=df)
sum.rel.lm <- summary(rel.lm)
print(sum.rel.lm)
df3$Si.Int.rel = round(predict(rel.lm, newdata = df2) + intercept, 5)

relSiPlot <- ggplot(df, aes(x=t.C.nm, y=Si.Int.rel)) +
    geom_line(color='red', size=1.25, aes(color="red"), data=df2) +
    annotate("text", label = 'LOESS',
           x = 25, y = 0.95,
           size = 5, colour = "red") +
    geom_line(color='blue', size=1.25, aes(color="blue"), data=df3) +
    annotate("text", label = 'lm',
           x = 25, y = 0.85,
           size = 5, colour = "blue") +
    geom_point(color='black', size=2) +
    xlab("C thickness [nm]") +
    ylab("Relative Si-K X-ray intensity ratio") +
    ggtitle("Monte Carlo Model of C on Si at 5 kV") +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14))#

print(relSiPlot)
fi <- '../mc3Scripts/dat/pdf/C-on-Si-5kV-10000-Traj-rel-si-predict.pdf'
ggsave(relSiPlot, file=fi, width=9.0, height=6.0, units="in", dpi=300 )

fi <- '../mc3Scripts/dat/png/C-on-Si-5kV-10000-Traj-rel-si-predict.png'
ggsave(relSiPlot, file=fi, width=9.0, height=6.0, units="in", dpi=300 )


print(head(df))

c.to.si.loess <- loess(c.to.si.mu ~ t.C.nm, span=0.75, data=df)

# let's save the fit...
fi <-'../mc3Scripts/dat/c-to-si-loess-5kV.RData'
save(c.to.si.loess, file=fi)

df$c.to.si.loess <-  predict(c.to.si.loess, data.frame(t.C.nm=df$t.C.nm))
# plot(df$t.C.nm, df$c.to.si.loess)


df2$c.to.si.mu = round(predict(c.to.si.loess, newdata = df2), 5)


ctosiInt <- ggplot(df, aes(x=t.C.nm, y=c.to.si.mu)) + 
    # coord_trans(y = "log") +
    geom_errorbar(aes(ymin=c.to.si.mu - 1.96*c.to.si.unc,
                      ymax=c.to.si.mu + 1.96*c.to.si.unc), width=.1) +
    geom_line(color='red', size=1.25, aes(color="red"), data=df2) +
    annotate("text", label = 'LOESS',
             x = 80, y = 1.25,
             size = 5, colour = "red") +
    geom_point(color='black', size=2) +
    xlab("C thickness [nm]") +
    ylab("C-K/Si-K X-ray intensity ratio") +
    ggtitle("Monte Carlo Model of C on Si at 5 kV") +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14))# or ,face="bold")) +



print(ctosiInt)

fi <- '../mc3Scripts/dat/pdf/C-on-Si-5kV-10000-Traj-predict.pdf'
ggsave(ctosiInt, file=fi, width=9.0, height=6.0, units="in", dpi=300 )

fi <- '../mc3Scripts/dat/png/C-on-Si-5kV-10000-Traj-predict.png'
ggsave(ctosiInt, file=fi, width=9.0, height=6.0, units="in", dpi=300 )

fi <-'../mc3Scripts/dat/csv/C-on-Si-5kV-10000-Traj-predict.csv'
write.csv(df2, file=fi, row.names=FALSE, quote=FALSE)
