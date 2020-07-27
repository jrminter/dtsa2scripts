# -*- coding: utf-8 -*-
#        1         2         3         4         5         6         7 |
# 3456789012345678901234567890123456789012345678901234567890123456789012
#
# simPdLineInCuMatrix.py
# jrm 2015-07-22 - Simulate the profile of a thin line embedded in a 
#                  block of another material. Nw main function is
#                  embedded in  dtsa2.jmMC3
#                  2:23:05.8 for 10000 electrons at 7 kV on crunch
import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3, FluorescenceXRayGeneration3, XRayTransport3", None)
import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.microanalysis.NISTMonte.Gen3 as nm3
import java.io as jio
import dtsa2.mcSimulate3 as mc3
import dtsa2.jmGen as jmg
import dtsa2.jmMC3 as jm3
import dtsa2 as dt2
import os
import shutil

# define the materials needed and sample size
# define the materials needed and sample size
 

umBlock  =    2.0  # size of the block in microns
nmLinWid =  100.0  # line width in nm
# define the electron optical parameters
e0       =   7.0 # kV
nTraj    = 100   # electrons
nPts     = 100   # number of points to compute
pc       =   1.5 # nA
lt       = 100.  # sec
dose     = lt*pc # nA*sec
imgSize =  512   # pixel size for images
imgSzUm  = 0.25
charF    = True  # include characteristic fluorescence
bremF    = True  # include continuum fluorescence 
poisN    = True  # include Poisson noise
vmrlEl   =  100  # number of el for VMRL

sc = 1.0e-06 # um to meters
umLine = nmLinWid * 1.0e-03

linMat = "Pd"
blkMat = "Cu"
lin = epq.Material(epq.Composition([epq.Element.Pd],[1.0],"Pd"), epq.ToSI.gPerCC(12.023))
blk = epq.Material(epq.Composition([epq.Element.Cu],[1.0],"Cu"), epq.ToSI.gPerCC(8.96))


det      = findDetector("Oxford p4 05eV 2K")
print(det)

homDir = os.environ['HOME']
relPrj = "/work/proj/QM15-04-07A-Ciminelli"
simDir = homDir + relPrj + "/dat/simDir"
csvDir = homDir + relPrj + "/dat/csv"
jmg.ensureDir(simDir)
jmg.ensureDir(csvDir)
# wd = homDir + relPrj + "/py/dtsa"
wd = homDir + relPrj + "/py/dtsa"
os.chdir(wd)
pyrDir = wd + "/simPdLineInCuMatrix Results"

#start clean
DataManager.clearSpectrumList()


# xrts=mc3.suggestTransitions("PdCu")

xrts = [epq.XRayTransition(epq.Element.Cu, epq.XRayTransition.LA1), epq.XRayTransition(epq.Element.Pd, epq.XRayTransition.LA1)]

xtraParams={}
xtraParams.update(mc3.configureXRayAccumulators(xrts, charAccum=charF, charFluorAccum=charF, bremFluorAccum=bremF))
xtraParams.update(mc3.configureOutput(simDir))
xtraParams.update(mc3.configureBeam(0.5*nmLinWid*1.0e-9, 0, -0.099, 1.0))
# xtraParams.update(mc3.configureGun(gun))
# mc3.useHeatMapPalette()

print(xtraParams)

# spc = jm3.lineInMatrix(lin, blk, nmLinWid, umBlock, det, e0, withPoisson=True, nTraj=nTraj, dose=dose, sf=charF, bf=bremF, xtraParams=xtraParams)
# spc = jm3.lineInMatrix(lin, blk, nmLinWid, umBlock, det, e0, poisN, nTraj, dose, charF, bremF, xtraParams)
#e0=20.0, dose=defaultDose, withPoisson=poisN, nTraj=defaultNumTraj, sf=defaultCharFluor, bf=defaultBremFluor, xtraParams=defaultXtraParams):
# spc = mc3.simulate(blk, det, e0, dose, poisN, nTraj, charF, bremF, xtraParams)
#
# Embedded Rectangle
#
# This works as expected.
# spc = mc3.embeddedRectangle(lin, [umLine*sc, umBlock*sc, umBlock*sc], blk, 0, det, e0, withPoisson=poisN, nTraj=nTraj, dose=dose, sf=charF, bf=bremF, xtraParams=xtraParams)
# 
# Multiblock
#
# spc = mc3.multiblock(blocks, det, e0, dose, True, nTraj, charF, bremF, xtraParams)
hdr ="x.um, kCuL, kPdL"
spc = jm3.simLineInMatrix3(lin, linMat, blk, blkMat, nmLinWid, umBlock, nPts, xrts, csvDir, hdr, det, e0, lt, pc, True, nTraj, True, True, 5, False, xtraParams={})
# spc.display()



shutil.rmtree(pyrDir)
print "Done!"