# -*- coding: utf-8 -*-
#         1         2         3         4         5         6         7 |
# 23456789012345678901234567890123456789012345678901234567890123456789012
#
# simCoatedOverBlock.py
# jrm 2015-10-02 - use mc3 to simulate a block of Ag on C covered with Al203
import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3, FluorescenceXRayGeneration3, XRayTransport3", None)
import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQLibrary.Detector as epd
import gov.nist.microanalysis.Utility as epu
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.microanalysis.NISTMonte.Gen3 as nm3
import gov.nist.microanalysis.EPQTools as et
import dtsa2.mcSimulate3 as mc3
import dtsa2.jmGen as jmg
import dtsa2.jmMC3 as jm3
import dtsa2 as dt2
import java.util as jutil
import java.io as jio
import java.nio.charset as cs
import os
import glob
import shutil

tFlm    =   0.20   # thickness um of film (Al2O3)
tBlk    =   5.0    # height um of block (Ag)
wBlk    =   5.0    # width um of block line
nTraj   =  250     # num Traj to run per pt 10000 for a long run
charF   =    True  # include characteristic fluorescence
bremF   =    True  # include continuum fluorescence 
pc      =    2.5   # nA
lt      =  100.0   # sec
e0      =   15.0   # keV
imgSize =  512     # pixel size for images
imgSzUm =   5.0    # image size in microns
vmrlEl  =  250     # number of el for VMRL
dose    = pc * lt  # nA sec
title   = "%g-um-Al2O3-on-%g-um-Ag-on-PET-%gkV" % (tFlm, tBlk, e0)

gitDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/mc3Scripts"
prjDir = gitDir + relPrj
datDir = prjDir + "/dat"
simDir = datDir + "/sim"
csvDir = datDir + "/csv"
msaDir = datDir + "/msa"

jmg.ensureDir(datDir)
jmg.ensureDir(csvDir)
jmg.ensureDir(simDir)
jmg.ensureDir(msaDir)

wd = gitDir + relPrj
os.chdir(wd)
pyrDir = wd + "/simCoatedOverBlock Results"

det  = findDetector("Oxford p4 05eV 2K")
# det  = findDetector("Si(Li)") # test with default
print(det)

pet   = epq.Material(epq.Composition([epq.Element.C,epq.Element.H,epq.Element.O],[0.62502,0.04196,0.069042]),epq.ToSI.gPerCC(1.37))
pet.setName("PET")
al2o3 = epq.Material(epq.Composition([epq.Element.Al,epq.Element.O],[0.5293,0.4707]),epq.ToSI.gPerCC(3.95))
al2o3.setName("Al2O3")
ag = epq.Material(epq.Composition([epq.Element.Ag],[1.0],"Ag"), epq.ToSI.gPerCC(11.9))
ag.setName("Ag")

# start clean
DataManager.clearSpectrumList()
# define the desired transitions (most intense for each line...)
#                     Ka1                  Ka1                   Ka1                      La1 
xrts = [transition("C K-L3"), transition("O K-L3"), transition("Al K-L3"), transition("Ag L3-M5")]
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

spc = jm3.coatedOverBlock(ag, tBlk*1.0e-6, wBlk*1.0e-6, al2o3, tFlm*1.0e-6, pet, det, e0=e0, withPoisson=True, nTraj=nTraj, dose=dose, sf=charF, bf=bremF, xtraParams=xtraParams)
display(spc)


shutil.rmtree(pyrDir)
print "Done!"