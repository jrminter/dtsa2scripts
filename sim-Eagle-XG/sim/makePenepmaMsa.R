library(rEDS)

homDir <- Sys.getenv('HOME')
homDir <- gsub('\\\\', '/', homDir)
relPath <- "/Documents/git/dtsa2Scripts/sim-Eagle-XG/sim"

inF <- sprintf("%s%s/%s", homDir, relPath,
               "20nm-C-on-EagleXG-Penepma/pe-spect-01.dat")

ouF <- sprintf("%s%s/%s", homDir, relPath,
               "pe-20nm-C-on-EagleXG-15kV.msa")

# inF <- "C:\\Users\\jrminter\\Desktop\\20nm-C-on-EagleXG\\pe-spect-01.dat"
# ouF <- "C:\\Users\\jrminter\\Desktop\\20nm-C-on-EagleXG\\pe-20nm-C-on-EagleXG-15kV.msa"

penepmaToMsa(inF, ouF, 15.0, 'PENEPMA 20 nm C on Eagle XG 15 kV')
