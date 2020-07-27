library(rEDS)

homDir <- Sys.getenv('HOME')
homDir <- gsub('\\\\', '/', homDir)
relPath <- "/Documents/git/dtsa2Scripts/sim-Al-vs-Br/sim"

inF <- sprintf("%s%s/%s", homDir, relPath,
               "Penepma-Al-5kV/pe-spect-01.dat")

ouF <- sprintf("%s%s/%s", homDir, relPath,
               "Penepma-Al-5kV.msa")

# inF <- "C:\\Users\\jrminter\\Desktop\\20nm-C-on-EagleXG\\pe-spect-01.dat"
# ouF <- "C:\\Users\\jrminter\\Desktop\\20nm-C-on-EagleXG\\pe-20nm-C-on-EagleXG-15kV.msa"

penepmaToMsa(inF, ouF, 5.0, 'PENEPMA Al 5 kV')

inF <- sprintf("%s%s/%s", homDir, relPath,
               "Penepma-KBr-5kV/pe-spect-01.dat")

ouF <- sprintf("%s%s/%s", homDir, relPath,
               "Penepma-KBr-5kV.msa")

penepmaToMsa(inF, ouF, 5.0, 'PENEPMA KBr 5 kV')
