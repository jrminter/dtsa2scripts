# -*- coding: utf-8 -*-
#         1         2         3         4         5         6         7 |
# 23456789012345678901234567890123456789012345678901234567890123456789012
#
# sim-C-on-Cu-w-Img.py
# jrm 2015-07-08 - use mc3 to simulate C on Cu
import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3, FluorescenceXRayGeneration3, XRayTransport3", None)
import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQLibrary.Detector as epd
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.microanalysis.NISTMonte.Gen3 as nm3
import gov.nist.microanalysis.EPQTools as et
import dtsa2.mcSimulate3 as mc3
import dtsa2.jmGen as jmg
import java.util as jutil
import java.io as jio
import java.nio.charset as cs
import os
import shutil

def ensureDir(d):
  """ensureDir(d)
  Check if the directory, d, exists, and if not create it."""
  if not os.path.exists(d):
    os.makedirs(d)

tNmC    =    25.0  # nm of C on Cu
nTraj   = 10000    # num Traj to run per pt 250 for a long run
charF   =    True  # include characteristic fluorescence
bremF   =    True  # include continuum fluorescence 
pc      =    2.5   # nA
lt      =  100.0   # sec
e0      =    7.0   # keV
imgSize =  512     # pixel size for images
imgSzUm =    1.0   # image size in microns
vmrlEl  =   40     # number of el for VMRL
dose = pc * lt # nA sec

gitDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/sim-C-on-Cu-w-Img"
datDir = gitDir + relPrj + "/msa"
csvDir = gitDir + relPrj + "/csv"
simDir = gitDir + relPrj + "/"
ensureDir(simDir)
ensureDir(datDir)

wd = gitDir + relPrj + "/py"
os.chdir(wd)
pyrDir = wd + "/sim-C-on-Cu-w-Img Results"

det  = findDetector("Oxford p4 05eV 2K")
print(det)


# start clean
DataManager.clearSpectrumList()

# create the materials
c  = epq.Material(epq.Composition([epq.Element.C], [1.0],"C"),  epq.ToSI.gPerCC(2.25))
cu = epq.Material(epq.Composition([epq.Element.Cu],[1.0],"Cu"), epq.ToSI.gPerCC(8.96))

# define the desired transitions
xrts=mc3.suggestTransitions("CCu")

# set up the extra parameters
xtraParams={}
xtraParams.update(mc3.configureXRayAccumulators(xrts, charAccum=charF, charFluorAccum=charF, bremFluorAccum=bremF))
# note that the image size on the specimen is in meters...
xtraParams.update(mc3.configureEmissionImages(xrts, imgSzUm*1.0e-6, imgSize))
xtraParams.update(mc3.configurePhiRhoZ(imgSzUm*1.0e-6))
xtraParams.update(mc3.configureTrajectoryImage(imgSzUm*1.0e-6, imgSize))
xtraParams.update(mc3.configureVRML(nElectrons = vmrlEl))
xtraParams.update(mc3.configureOutput(simDir))

print(xtraParams)

# treat the substrate as 10 um C
sLay = [[c, tNmC/1.0e9], [cu, 10.0/1.0e6]]
sSpc = mc3.multiFilm(sLay, det, e0=e0, withPoisson=True, nTraj=nTraj, dose=dose, sf=charF, bf=bremF, xtraParams=xtraParams)
display(sSpc)


# clean up cruft
shutil.rmtree(pyrDir)
print "Done!"

