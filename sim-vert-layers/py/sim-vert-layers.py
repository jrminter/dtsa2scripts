# -*- coding: utf-8 -*-
#         1         2         3         4         5         6         7 |
# 23456789012345678901234567890123456789012345678901234567890123456789012
#
# sim-vert-layers.py
# jrm 2015-07-09 - use mc3 to simulate vertical layers
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
import glob
import shutil

widNm   =   25.0   # width of layers [nm]
nTraj   =  250     # num Traj to run per pt 10000 for a long run
charF   =    True  # include characteristic fluorescence
bremF   =    True  # include continuum fluorescence 
pc      =    2.5   # nA
lt      =  100.0   # sec
e0      =    7.0   # keV
imgSize =  512     # pixel size for images
imgSzUm =   5.0   # image size in microns
vmrlEl  =   40     # number of el for VMRL
dose    = pc * lt # nA sec


gitDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/sim-vert-layers"
simDir = gitDir + relPrj + "/"
jmg.ensureDir(simDir)

wd = gitDir + relPrj + "/py"
os.chdir(wd)
pyrDir = wd + "/sim-vert-layers Results"

det  = findDetector("Oxford p4 05eV 2K")
print(det)


# start clean
DataManager.clearSpectrumList()

pet = epq.Material(epq.Composition([epq.Element.C,epq.Element.H,epq.Element.O],[0.62502,0.04196,0.069042]),epq.ToSI.gPerCC(1.37))
cu = epq.Material(epq.Composition([epq.Element.Cu],[1.0],"Cu"), epq.ToSI.gPerCC(8.96))
pd = epq.Material(epq.Composition([epq.Element.Pd],[1.0],"Pd"), epq.ToSI.gPerCC(11.9))

# define the desired transitions
# xrts = [epq.XRayTransitionSet(epq.Element.Pd, epq.XRayTransitionSet.L_FAMILY), epq.XRayTransitionSet(epq.Element.Cu, epq.XRayTransitionSet.L_FAMILY)]
xrts=mc3.suggestTransitions("CuPd")
print(xrts)
for tr in xrts:
  print(tr.getSiegbahnName())
for tr in xrts:
  print(tr.getIUPACName())
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

spc = mc3.verticalLayers(widNm*1.0e-9, [pd,cu], det, e0=e0, withPoisson=True, nTraj=nTraj, dose=dose, sf=charF, bf=bremF, xtraParams=xtraParams)
display(spc)

shutil.rmtree(pyrDir)
print "Done!"