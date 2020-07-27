# -*- coding: utf-8 -*-
#         1         2         3         4         5         6         7 |
# 23456789012345678901234567890123456789012345678901234567890123456789012
#
# sim-SiC.py
# jrm 2019-12-14 - use mc3 to simulate SiC for test Pouchou


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

gitDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/sim-Pouchou-SiC"
datDir = gitDir + relPrj + "/msa"
csvDir = gitDir + relPrj + "/csv"
ensureDir(csvDir)
wd = gitDir + relPrj + "/py"
os.chdir(wd)
pyrDir = wd + "/sim-Pouchou-SiC"
ensureDir(datDir)
det  = findDetector("Oxford p4 05eV 2K")
print(det)

nTraj    =  1000   # num Traj to run per pt 10000 for a long run
charF    =  True   # include characteristic fluorescence
bremF    =  True   # include continuum fluorescence 
pc       =     2.5  # nA
lt       =   100.0  # sec
e0       =     3.0  # keV

dose = pc * lt # nA sec
# start clean
DataManager.clearSpectrumList()

# create the materials
si       = material("Si", density=2.32)
sic      = material("SiC", density=3.21)
y3fe5o12 = material("Y3Fe5O12", density=5.11)

trs = [epq.XRayTransitionSet(epq.Element.C, epq.XRayTransitionSet.K_ALPHA),
       epq.XRayTransitionSet(epq.Element.O, epq.XRayTransitionSet.K_ALPHA),
       epq.XRayTransitionSet(epq.Element.Si, epq.XRayTransitionSet.K_ALPHA)]

# start with bulk standards
#
# C

siSpc = mc3.simulate(si, det, e0=e0, dose=dose, withPoisson=True, nTraj=nTraj, sf=charF, bf=bremF, xtraParams={})
sp=siSpc.getProperties()
sp.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName,"Si-std")
sp.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
sp.setNumericProperty(epq.SpectrumProperties.FaradayBegin, dose)
sp.setNumericProperty(epq.SpectrumProperties.FaradayEnd, dose)
sp.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
sp.setCompositionProperty(epq.SpectrumProperties.StandardComposition, epq.Composition(epq.Element.Si))
display(siSpc)
siStd  = {"El":element("Si"),  "Spc":siSpc}

sic = epq.Material(epq.Composition([epq.Element.Si,epq.Element.C],
                                   [0.70045,0.29955]),epq.ToSI.gPerCC(3.21))
sic.setName("SiC")

sicSpc = mc3.simulate(sic, det, e0=e0, dose=dose, withPoisson=True, nTraj=nTraj, sf=charF, bf=bremF, xtraParams={})
sp=sicSpc.getProperties()
sp.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName,"SiC-std")
sp.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
sp.setNumericProperty(epq.SpectrumProperties.FaradayBegin, dose)
sp.setNumericProperty(epq.SpectrumProperties.FaradayEnd, dose)
sp.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
sp.setCompositionProperty(epq.SpectrumProperties.StandardComposition, sic)
display(sicSpc)
# sicStd  = {"El":e,  "Spc":siSpc}

y3fe5o12Spc = mc3.simulate(y3fe5o12, det, e0=e0, dose=dose, withPoisson=True, nTraj=nTraj, sf=charF, bf=bremF, xtraParams={})
sp=y3fe5o12Spc.getProperties()
sp.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName,"Y3Fe5O12-std")
sp.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
sp.setNumericProperty(epq.SpectrumProperties.FaradayBegin, dose)
sp.setNumericProperty(epq.SpectrumProperties.FaradayEnd, dose)
sp.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
sp.setCompositionProperty(epq.SpectrumProperties.StandardComposition, y3fe5o12)
display(y3fe5o12Spc)
# y3fe5o12Std  = {"El":element("Si"),  "Spc":siSpc}


"""

stds = [cuStd, auStd]

# list of thicknesses
lNm=[]
# lists of Au K-ratios
lKAuL=[]
# lKAuM=[]
# lists of Cu K-ratios
# lKCuL=[]
lKCuK=[]

for i in range(nSteps):
    tNmAu = tNmStep*float(i+1)
    print(tNmAu)
    lNm.append(tNmAu)
    sLay = [[au, tNmAu/1.0e9], [cu, 100.0]]
    sSpc = mc3.multiFilm(sLay, det, e0=e0, withPoisson=True, nTraj=nTraj, dose=dose, sf=charF, bf=bremF, xtraParams={})
    sp=sSpc.getProperties()
    sp.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName,"%g-nm-Au-on-Cu"% tNmAu)
    sp.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
    sp.setNumericProperty(epq.SpectrumProperties.FaradayBegin, dose)
    sp.setNumericProperty(epq.SpectrumProperties.FaradayEnd, dose)
    display(sSpc)
    a = jmg.compKRs(sSpc, stds, trs, det, e0)
    print(a[0])
    print(a[1])
    lKAuL.append(a[0])
    lKCuK.append(a[1])
    print (i+1)


# prepare the output file
csvPath = csvDir + "/sim-Au-on-Cu-%g-kV-%g-steps-%g-nm.csv" % (e0, nSteps, tNmStep)
f=open(csvPath, 'w')
strLine = 'tNm, kAuLaMu, kAuLaStd, kCuKaMu, kCuKaStd\n'
f.write(strLine)
l = len(lNm)
for i in range(l):
    strLine = '%.1f, %.5f, %.5f, %.5f, %.5f\n' % (lNm[i], lKAuL[i][0], lKAuL[i][1], lKCuK[i][0], lKCuK[i][1])
    f.write(strLine)

f.close()

"""

# clean up cruft
# shutil.rmtree(pyrDir)
print "Done!"




