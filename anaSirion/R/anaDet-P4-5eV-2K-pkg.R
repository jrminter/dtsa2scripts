# anaDet-P4-5eV-2K.R
rm(list=ls())

library(ggplot2)
library(rEDS)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

gitDir <- Sys.getenv("GIT_HOME")
rDir <- paste0(gitDir, '/dtsa2Scripts/anaSirion/R/')
setwd(rDir)


csvFil <- '../dat/Oxford-P4-05eV-2K-15kV-Cu.csv'
res <- anaDetectorFile(csvFil, "Oxford-P4-05eV-2K", labOffCW=3,
                       labOffZO=3, labOffMn=3, labOffCu=3)
# print(res) #['detCW'])

library(grid)
library(gridExtra)
p1 <- res['detCW']
p2 <- res['detZO']
p3 <- res['detMn']
p4 <- res['cuLaPk']

print(class(p1))

grid.arrange(p1, p2, p3, p4, ncol = 2, main = "Main title")
# multiplot(res['detCW'], res['detZO'], res['detMn'], res['cuLaPk'], cols=2)
