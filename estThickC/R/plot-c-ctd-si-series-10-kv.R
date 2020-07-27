# plot-c-ctd-si-series.R

rm(list=ls())
library(ggplot2)
library(grid)
library(gridExtra)

gitDir <- Sys.getenv("GIT_HOME")
str.wd <- paste0(gitDir, '/dtsa2Scripts/estThickC/R/')
setwd(str.wd)

fi <-'../dat/csv/C-ctd-Si-10-kV-10000-Traj.csv'

df <- read.csv(fi, header=TRUE, as.is=TRUE)
print(head(df))

cInt <- ggplot(df, aes(x=t.C.nm, y=C.Int.mu)) + 
  geom_errorbar(aes(ymin=C.Int.mu - 1.96*C.Int.unc,
                    ymax=C.Int.mu + 1.96*C.Int.unc), width=.1) +
  # geom_line(color='blue') +
  geom_point(color='blue') + xlab("C thicknesss [nm]") +
  ylab("C-K intensity") #+
  #annotate("text", x = -3.75, y = 1750, label = "shell") +
  #annotate("text", x = 0, y = 1750, label = "core") +
  #annotate("text", x = 3.75, y = 1750, label = "shell")

print(cInt)


siInt <- ggplot(df, aes(x=t.C.nm, y=Si.Int.mu)) + 
  geom_errorbar(aes(ymin=Si.Int.mu - 1.96*Si.Int.unc,
                    ymax=Si.Int.mu + 1.96*Si.Int.unc), width=.1) +
  # geom_line(color='blue') +
  geom_point(color='blue') + xlab("C thicknesss [nm]") +
  ylab("Si-K intensity") #+
#annotate("text", x = -3.75, y = 1750, label = "shell") +
#annotate("text", x = 0, y = 1750, label = "core") +
#annotate("text", x = 3.75, y = 1750, label = "shell")

print(siInt)

df$c.to.si.mu <- df$C.Int.mu / df$Si.Int.mu
df$c.to.si.unc <- sqrt((df$C.Int.unc/df$C.Int.mu)**2 + 
                        (df$Si.Int.unc/df$Si.Int.mu)**2)



c.to.si.loess.10kV <- loess(c.to.si.mu~t.C.nm, df)
print(summary(c.to.si.loess.10kV))
fi <- '../dat/rda/c-to-si-loess-10kV.rda'
save(c.to.si.loess.10kV, file=fi)


df2 <- df[-1,]

exp.model <-lm(log(c.to.si.mu)~ t.C.nm, data=df2) 
lin.model <- lm(c.to.si.mu~ 0+ t.C.nm, data=df2)
print(summary(lin.model))



t.C.nm <- seq(from=0.5, to=100, by=0.5)
df3 <- data.frame(t.C.nm=t.C.nm)
df3$c.to.si.mu = exp(predict(exp.model, newdat=df3))

poly.model <- lm(c.to.si.mu ~ 0 + I(t.C.nm^2) + I(t.C.nm^3) +
                 I(t.C.nm^4) + I(t.C.nm^5) + I(t.C.nm^6),
                 data=df2)

df4 <- data.frame(t.C.nm=t.C.nm)
df4$c.to.si.mu = predict(poly.model, newdat=df4)

ss6 <- smooth.spline(df4$t.C.nm, df4$c.to.si.mu, df = 10)

ctosiInt <- ggplot(df2, aes(x=t.C.nm, y=c.to.si.mu)) + 
    # coord_trans(y = "log") +
    geom_errorbar(aes(ymin=c.to.si.mu - 1.96*c.to.si.unc,
                      ymax=c.to.si.mu + 1.96*c.to.si.unc), width=.1) +
    # geom_line(color='red', size=1.25, aes(color="red"), data=df4) +
    # geom_smooth() +
    stat_smooth(method = "loess") +
    # stat_smooth( method = "lm", formula = y ~ 0 + x) +
    annotate("text", label = 'loess', # loess smooth',
             x = 70, y = 0.09,
             size = 5, colour = "blue") +
    geom_point(color='black', size=2) +
    xlab("C thickness [nm]") +
    ylab("C-K/Si-K X-ray intensity ratio") +
    ggtitle("Monte Carlo Model of C-coatings on Si at 10 kV") +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          plot.title = element_text(hjust = 0.5)) # center the title

# df4$c.to.si.mu = predict(loess(c.to.si.mu~t.C.nm, df4))
df4$c.to.si.mu = predict(lm(c.to.si.mu~0 + t.C.nm, df4))
plot(c.to.si.mu~t.C.nm, data=df4)


print(ctosiInt)

fi <- '../pdf/c-ctd-si-series-plt-10kV.pdf'
ggsave(ctosiInt, file=fi, width=9.0, height=6.0, units="in", dpi=300)

fi <- '../png/c-ctd-si-series-plt-10kV.png'
ggsave(ctosiInt, file=fi, width=9.0, height=6.0, units="in", dpi=113.7778)

fi <- '../dat/csv/c-ctd-si-loess-pred-10kV.csv'
write.csv(df4, file=fi, row.names=FALSE)
