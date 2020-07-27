# -*- coding: utf-8 -*-
#         1         2         3         4         5         6         7 |
# 23456789012345678901234567890123456789012345678901234567890123456789012
#
# sim-C-on-Cu.py
# jrm 2015-07-04 - use mc3 to simulate C on Cu
#  7 kV Elapse: 0:08:57.1 on jrmFastMac
# 10 kV         0:13:22.3 on jrmFastMac
# 15 kV         0:21:22.0 on jrmFastMac
# 20 kV         0:28:51.5 on jrmFastMac
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
relPrj = "/dtsa2Scripts/sim-C-on-Cu"
datDir = gitDir + relPrj + "/msa"
csvDir = gitDir + relPrj + "/csv"
ensureDir(csvDir)
wd = gitDir + relPrj + "/py"
os.chdir(wd)
pyrDir = wd + "/sim-C-on-Cu Results"
ensureDir(datDir)
det  = findDetector("Oxford p4 05eV 2K")
print(det)

tNmC    =   500.0  # nm C to sim
nTraj   =   100    # num Traj to run per pt 250 for a long run
charF   =   True   # include characteristic fluorescence
bremF   =   True   # include continuum fluorescence 
pc      =   2.5    # nA
lt      = 100.0    # sec
e0      =   7.0    # keV

dose = pc * lt # nA sec
# start clean
DataManager.clearSpectrumList()

# create the materials
c  = material("C", density=2.25)
cu = material("Cu", density=8.96)

trs = [epq.XRayTransitionSet(epq.Element.C, epq.XRayTransitionSet.K_FAMILY), 
epq.XRayTransitionSet(epq.Element.Cu, epq.XRayTransitionSet.L_FAMILY)]

# start with bulk standards

cSpc = mc3.simulate(c, det, e0=e0, dose=dose, withPoisson=True, nTraj=nTraj, sf=charF, bf=bremF, xtraParams={})
sp=cSpc.getProperties()
sp.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName,"C-std")
sp.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
sp.setNumericProperty(epq.SpectrumProperties.FaradayBegin, dose)
sp.setNumericProperty(epq.SpectrumProperties.FaradayEnd, dose)
sp.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
sp.setCompositionProperty(epq.SpectrumProperties.StandardComposition, epq.Composition(epq.Element.C))
# display(cSpc)
cStd  = {"El":element("C"),  "Spc":cSpc}

cuSpc = mc3.simulate(cu, det, e0=e0, dose=dose, withPoisson=True, nTraj=nTraj, sf=charF, bf=bremF, xtraParams={})
sp=cuSpc.getProperties()
sp.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName,"Cu-std")
sp.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
sp.setNumericProperty(epq.SpectrumProperties.FaradayBegin, dose)
sp.setNumericProperty(epq.SpectrumProperties.FaradayEnd, dose)
sp.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
sp.setCompositionProperty(epq.SpectrumProperties.StandardComposition, epq.Composition(epq.Element.Cu))
# display(cuSpc)
cuStd = {"El":element("Cu"), "Spc":cuSpc}
stds = [cStd, cuStd]

lNm=[]
lKCK=[]
lKCuL=[]

for i in range(500):
    tNmC = float(i+1)
    lNm.append(tNmC)
    sLay = [[c, tNmC/1.0e9], [cu, 100.0]]
    sSpc = mc3.multiFilm(sLay, det, e0=e0, withPoisson=True, nTraj=nTraj, dose=dose, sf=charF, bf=bremF, xtraParams={})
    sp=sSpc.getProperties()
    sp.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName,"%g-nm-C-on-Cu"% tNmC)
    sp.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
    sp.setNumericProperty(epq.SpectrumProperties.FaradayBegin, dose)
    sp.setNumericProperty(epq.SpectrumProperties.FaradayEnd, dose)
    # display(sSpc)
    a = jmg.compKRs(sSpc, stds, trs, det, e0)
    lKCK.append(a[0])
    lKCuL.append(a[1])
    print (i+1)

# prepare the output file
csvPath = csvDir + "/sim-C-on-Cu-%gkV.csv" % e0
f=open(csvPath, 'w')
strLine = 'tNm, kCKa, kCuLa\n'
f.write(strLine)
l = len(lNm)
for i in range(l):
    strLine = '%.1f, %.5f, %.5f\n' % (lNm[i], lKCK[i], lKCuL[i])
    f.write(strLine)

f.close()



# clean up cruft
shutil.rmtree(pyrDir)
print "Done!"




