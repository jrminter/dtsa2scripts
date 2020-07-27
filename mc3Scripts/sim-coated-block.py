# -*- coding: utf-8 -*-
#         1         2         3         4         5         6         7 |
# 23456789012345678901234567890123456789012345678901234567890123456789012
#
# sim-coated-block.py
# jrm 2015-09-22 - use mc3 to simulate a block of Ag on C covered with Al203
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

tFlm   =   0.10    # thickness um of film (Al2O3)
tBlk   =   0.50    # height um of block (Ag)
wBlk    =  0.50    # width um of block line
nTraj   =  200     # num Traj to run per pt 10000 for a long run
charF   =    True  # include characteristic fluorescence
bremF   =    True  # include continuum fluorescence 
pc      =    2.5   # nA
lt      =  100.0   # sec
e0      =   15.0   # keV
imgSize =  512     # pixel size for images
imgSzUm =   5.0    # image size in microns
vmrlEl  =   40     # number of el for VMRL
dose    = pc * lt  # nA sec
title   = "%g-um-Al2O3-on-%g-um-Ag-on-PET-%gkV" % (tFlm, tBlk, e0)

gitHom = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/mc3Scripts"
prjDir = gitHom + relPrj
datDir = prjDir + "/dat"
csvDir = datDir + "/csv"
simDir = datDir + "/sim"

jmg.ensureDir(datDir)
jmg.ensureDir(csvDir)
jmg.ensureDir(simDir)


os.chdir(prjDir)
pyrDir = prjDir + "/sim-coated-block Results"

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

def coatedOverBlock(mat, height, width, coating, thickness, substrate, det, e0=20.0, withPoisson=True, nTraj=defaultNumTraj, dose=defaultDose, sf=defaultCharFluor, bf=defaultBremFluor, xtraParams=defaultXtraParams):
   """coatedOverBlock(mat, height, width, coating, thickness, substrate, det, e0=20.0, withPoisson=True, nTraj=defaultNumTraj, dose=defaultDose, sf=defaultCharFluor, bf=defaultBremFluor, substrate=None, xtraParams={})
   Monte Carlo simulate a spectrum from a block shaped particle of the specified material (mat) and height (z in m) and width (x and y in m). \
   The block and subtrate is coated in a material 'coating' of the specified thickness which fully encapsulates the particle and covers the substrate too."""
   def buildBlock(monte, chamber, origin, buildParams):
      height = buildParams["Height"]
      width = buildParams["Width"]
      subMat = buildParams["Substrate"]
      mat = buildParams["Material"]
      coating = buildParams["Coating"]
      thickness = buildParams["Thickness"]
      
      ccDims = [width + 2.0 * thickness, width + 2.0 * thickness, height + thickness]
      ccPt   = epu.Math2.plus(origin, [0.0, 0.0, 0.5*(height + thickness)])
      coatedCube = nm.MultiPlaneShape.createBlock(ccDims, ccPt, 0.0, 0.0, 0.0)
      sr1 = monte.addSubRegion(chamber, coating, coatedCube)

      cDims = [width, width, height]
      cPt   = epu.Math2.plus(origin, [0.0, 0.0, thickness+0.5*height])
      cube = nm.MultiPlaneShape.createBlock(cDims, cPt, 0.0, 0.0, 0.0)

      monte.addSubRegion(sr1, mat, cube)

      sideSlabWidth = 2.5*1.0e-6 - (thickness + 0.5*width)
      sideSlabDims = [sideSlabWidth, 5.*1.0e-6, thickness]
      leftSidePt = epu.Math2.plus(origin, [0.5*(width+sideSlabWidth), 0.0, thickness+height])
      leftSide = nm.MultiPlaneShape.createBlock(sideSlabDims, leftSidePt, 0.0, 0.0, 0.0)
      monte.addSubRegion(chamber, coating, leftSide)
      rightSidePt = epu.Math2.plus(origin, [-0.5*(width+sideSlabWidth), 0.0, thickness+height])
      rightSide = nm.MultiPlaneShape.createBlock(sideSlabDims, rightSidePt, 0.0, 0.0, 0.0)
      monte.addSubRegion(chamber, coating, rightSide)

      subNorm = [0.0, 0.0, -1.0]
      subPt   = epu.Math2.plus(origin, [0.0, 0.0, height + thickness])
      monte.addSubRegion(chamber, subMat, nm.MultiPlaneShape.createSubstrate(subNorm, subPt))

   s1 = u"MC simulation of a [%0.2f,%0.2f,%0.2f] micron block of %s%s" % (width * 1.0e6, width * 1.0e6, height * 1.0e6, mat, (" on %s" % substrate if substrate else ""))
   s2 = u" coated with %0.2f microns of %s at %0.1f keV%s%s" % (thickness* 1.0e6, coating, e0, (" + CSF" if sf else ""), (" + BSF" if bf else ""))
   tmp = s1 + s2

   # tmp = u"MC simulation of a [%0.2f,%0.2f,%0.2f] micron block of %s%s coated with %0.2f microns of %s at %0.1f keV%s%s" % (width * 1.0e6, width * 1.0e6, height * 1.0e6, mat, (" on %s" % substrate if substrate else ""), coating, e0, (" + CSF" if sf else ""), (" + BSF" if bf else ""))
   params = {"Substrate": substrate, "Width" : width, "Height" : height, "Material" : mat, "Coating" : coating, "Thickness" : thickness}
   return mc3.base(det, e0, withPoisson, nTraj, dose, sf, bf, tmp, buildBlock, params, xtraParams)


pet   = epq.Material(epq.Composition([epq.Element.C,epq.Element.H,epq.Element.O],[0.62502,0.04196,0.069042]),epq.ToSI.gPerCC(1.37))
pet.setName("PET")
al2o3 = epq.Material(epq.Composition([epq.Element.Al,epq.Element.O],[0.5293,0.4707]),epq.ToSI.gPerCC(3.95))
al2o3.setName("Al2O3")
ag = epq.Material(epq.Composition([epq.Element.Ag],[1.0],"Ag"), epq.ToSI.gPerCC(11.9))
ag.setName("Ag")

# start clean
DataManager.clearSpectrumList()
# define the desired transitions
xrts=mc3.suggestTransitions("COAlAg")
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

spc = coatedOverBlock(ag, tBlk*1.0e-6, wBlk*1.0e-6, al2o3, tFlm*1.0e-6, pet, det, e0=e0, withPoisson=True, nTraj=nTraj, dose=dose, sf=charF, bf=bremF, xtraParams=xtraParams)
display(spc)


shutil.rmtree(pyrDir)
print "Done!"