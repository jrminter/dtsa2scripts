# -*- coding: utf-8 -*-
#         1         2         3         4         5         6         7 |
# 23456789012345678901234567890123456789012345678901234567890123456789012
#
# sim-Au-on-Cu.py
# jrm 2019-11-14 - use mc3 to simulate Au film on Cu
#  15 kV Elapse: 0:08:57.1 on jrmFastMac
# jrm 2019-12-03 - trying 15 steps of 20 nm with 10,000 trajectories on
#                  jrmSimulationPC
#                  Starting at 9:07 p.m Elapse: 0:51:32.7
# jrm 2019-12-04 - trying 30 steps of 10 nm with 10,000 trajectories on
#                  jrmSimulationPC
#                  Elapse: 1:34:27.4


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
relPrj = "/dtsa2Scripts/sim-Au-on-Cu"
datDir = gitDir + relPrj + "/msa"
csvDir = gitDir + relPrj + "/csv"
ensureDir(csvDir)
wd = gitDir + relPrj + "/py"
os.chdir(wd)
pyrDir = wd + "/sim-Au-on-Cu Results"
ensureDir(datDir)
det  = findDetector("Oxford p4 05eV 4K")
print(det)

tNmStep  =    10.0  # nm Au step size
nSteps   =    30    # Num
nTraj    = 100    # num Traj to run per pt 10000 for a long run
charF    =   True   # include characteristic fluorescence
bremF    =   True   # include continuum fluorescence 
pc       =     2.5  # nA
lt       =   100.0  # sec
e0       =    15.0  # keV

dose = pc * lt # nA sec
# start clean
DataManager.clearSpectrumList()

# create the materials
ag = material("Ag", density=10.5)


trs = [epq.XRayTransitionSet(epq.Element.Au, epq.XRayTransitionSet.L_ALPHA),
epq.XRayTransitionSet(epq.Element.Cu, epq.XRayTransitionSet.K_ALPHA)]

# start with bulk standards
#
# Au

auSpc = mc3.simulate(au, det, e0=e0, dose=dose, withPoisson=True, nTraj=nTraj, sf=charF, bf=bremF, xtraParams={})
sp=auSpc.getProperties()
sp.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName,"Au-std")
sp.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
sp.setNumericProperty(epq.SpectrumProperties.FaradayBegin, dose)
sp.setNumericProperty(epq.SpectrumProperties.FaradayEnd, dose)
sp.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
sp.setCompositionProperty(epq.SpectrumProperties.StandardComposition, epq.Composition(epq.Element.Au))
display(auSpc)
auStd  = {"El":element("Au"),  "Spc":auSpc}

cuSpc = mc3.simulate(cu, det, e0=e0, dose=dose, withPoisson=True, nTraj=nTraj, sf=charF, bf=bremF, xtraParams={})
sp=cuSpc.getProperties()
sp.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName,"Cu-std")
sp.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
sp.setNumericProperty(epq.SpectrumProperties.FaradayBegin, dose)
sp.setNumericProperty(epq.SpectrumProperties.FaradayEnd, dose)
sp.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
sp.setCompositionProperty(epq.SpectrumProperties.StandardComposition, epq.Composition(epq.Element.Cu))
display(cuSpc)
cuStd = {"El":element("Cu"), "Spc":cuSpc}


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


# clean up cruft
shutil.rmtree(pyrDir)
print "Done!"




