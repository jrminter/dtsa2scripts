# -*- coding: utf-8 -*-
#         1         2         3         4         5         6         7 |
# 23456789012345678901234567890123456789012345678901234567890123456789012
#
# sim-bilayer-on-sub.py
# jrm 2015-07-13 - use mc3 to simulate a bilayer line
import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3, FluorescenceXRayGeneration3, XRayTransport3", None)
import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQLibrary.Detector as epd
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.microanalysis.NISTMonte.Gen3 as nm3
import gov.nist.microanalysis.EPQTools as et
import dtsa2.mcSimulate3 as mc3
import dtsa2.jmGen as jmg
import dtsa2.jmMC3 as jm3
import java.util as jutil
import java.io as jio
import java.nio.charset as cs
import os
import glob
import shutil

tNmCu   =  300.0   # nm of Cu top layer
tNmAg   =   40.0   # nm of Ag top layer
wUm     =    5.0   # width microns of line
lUm     =    5.0   # length microns of line
nTraj   =  200     # num Traj to run per pt 10000 for a long run
charF   =    True  # include characteristic fluorescence
bremF   =    True  # include continuum fluorescence 
pc      =    2.5   # nA
lt      =  100.0   # sec
e0      =   30.0   # keV
imgSize =  512     # pixel size for images
imgSzUm =   5.0   # image size in microns
vmrlEl  =   40     # number of el for VMRL
dose    = pc * lt # nA sec
title   = "%gnm-Cu-on-%gAg-on-PET-%gkV" % (tNmCu, tNmAg, e0)

gitDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/sim-bilayer-on-sub"
datDir = gitDir + relPrj + "/msa"
csvDir = gitDir + relPrj + "/csv"
simDir = gitDir + relPrj + "/"
jmg.ensureDir(simDir)
jmg.ensureDir(datDir)

wd = gitDir + relPrj + "/py"
os.chdir(wd)
pyrDir = wd + "/sim-bilayer-on-sub Results"

det  = findDetector("Oxford p4 05eV 2K")
print(det)

if 'defaultXtraParams' not in globals():
   defaultXtraParams = {}
if 'defaultBremFluor' not in globals():
   defaultBremFluor = False
if 'defaultCharFluor' not in globals():
   defaultCharFluor = False
if 'defaultNumTraj' not in globals():
   defaultNumTraj = 1000
if 'defaultDose' not in globals():
   defaultDose = 120.0


# start clean
DataManager.clearSpectrumList()


pet = epq.Material(epq.Composition([epq.Element.C,epq.Element.H,epq.Element.O],[0.62502,0.04196,0.069042]),epq.ToSI.gPerCC(1.37))
cu = epq.Material(epq.Composition([epq.Element.Cu],[1.0],"Cu"), epq.ToSI.gPerCC(8.96))
ag = epq.Material(epq.Composition([epq.Element.Ag],[1.0],"Ag"), epq.ToSI.gPerCC(11.9))

# define the desired transitions
xrts=mc3.suggestTransitions("COCuAg")
# print(xrts)
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

spc = jm3.biLayerLineOnSubstrate(cu, ag, pet, tNmCu*1.0e-9, tNmAg*1.0e-9, wUm*1.0e-06, lUm*1.0e-06, det, title, e0=e0, withPoisson=True, nTraj=nTraj, dose=dose, sf=charF, bf=bremF, xtraParams=xtraParams)
display(spc)

name = "biLayerLineOnSub-%s" % title
strFi = simDir + "/" + name +"/" + name + ".msa"

spc.save(strFi)

shutil.rmtree(pyrDir)
print "Done!"