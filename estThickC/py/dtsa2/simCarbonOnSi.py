# -*- coding: utf-8 -*-
#         1         2         3         4         5         6         7 |
# 23456789012345678901234567890123456789012345678901234567890123456789012
#
# simCarbonOnSi.py
# JRM 2016-08-06
# revised 2017-07-17
# 
# Use MC3 coated substrate
# Elapse: 0:31:53.3 on crunch
#         1:01:01.4 on ROCPW7C5C42 at 10 kV
#         1:50:33.1 on ROCPW7C5C42 at 20 kV
#         1:47:58.8 on crunch at 7 kV for 2 nm steps to 100 nm
import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3,FluorescenceXRayGeneration3, XRayTransport3", None)

import os
import glob
import shutil
import time
import math
import csv

import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.microanalysis.NISTMonte.Gen3 as nm3
import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQLibrary.Detector as epd
import gov.nist.microanalysis.Utility as epu
import gov.nist.microanalysis.EPQTools as ept


import dtsa2 as dt2
import dtsa2.mcSimulate3 as mc3

import dtsa2.jmGen as jmg

start    = time.time()

e0       = 7.0    # beam voltage(kV)
nTraj    = 10000  # number of trajectories to compute
pc       =  1.0   # probe current (nA)
lt       = 100    # spectrum live time (sec)
detNam   = "Oxford p4 05eV 2K"  #  detector name
imgSzPx  = 512
imgSzUm  = 2.0
beamSzNm = 1.0
tNmMaxC  = 100
tNmStepC = 2


homDir  = os.environ['GIT_HOME']
relDir  = "/dtsa2Scripts/estThickC"

inDir   = homDir + relDir + "/py/dtsa2"
outDir  = homDir + relDir + "/dat/csv"

jmg.ensureDir(outDir)


rptDir  = inDir + '/simCarbonOnSi Results/'

def anaConSi(spc, det, digits=2, display=False):
    """anaConSi

    Strip continuum and analyze a spectrum for C, on CuL

    Parameters
    ----------
    spc: A DTSA-II scriptable spectrum
        The spectrum to process
    det: The detector

    Returns
    -------
    A dictionary of tuples for C and Cu where each tuple is the
    peak integral (in counts/nA-sec)  and the uncertainty.
    """
    sp = spc.getProperties()
    pc = sp.getNumericProperty(epq.SpectrumProperties.FaradayBegin)
    lt = sp.getNumericProperty(epq.SpectrumProperties.LiveTime)
    dose = pc*lt
    sp.setDetector(det)
    bks = jmg.StripBackground(spc, det)
    if display:
        nam = spc.getProperties().getTextProperty(epq.SpectrumProperties.SpectrumDisplayName) + '-bks'
        bks.rename(nam)
        bks.display()
    cInt  = jmg.compPeakIntegral(bks, 0.282, 156, 1)
    siInt = jmg.compPeakIntegral(bks, 1.746, 240, 1)

    out = { "C" : (round(cInt[0]/dose, digits), round(cInt[1]/dose, digits)),
            "Si": (round(siInt[0]/dose, digits), round(siInt[1]/dose, digits)) }

    return out


c = material("C", density=2.1)
si = material("Si", density=2.3296)

det = findDetector(detNam)
dose = lt*pc

DataManager.clearSpectrumList()

lT       = [] # will have t=0...
lPICmu   = [] # an array for mean C peak integral
lPICuc   = [] # an array for C peak integral uncertainty
lPISimu  = [] # an array for mean Cu peak integral
lPISiuc  = [] # an array for Cu peak integral uncertainty

# empty xtra-param
xp = {}

# empty x-ray transition set
xrts = []
# set a counter
iCount = 0

# get and add the C-K family
cTS  = epq.XRayTransitionSet(epq.Element.C,
                             epq.XRayTransitionSet.K_FAMILY)
cTrs = cTS.getTransitions()
for tr in cTrs:
    xrts.append(tr)

# get and add the Cu-L family
siTS  = epq.XRayTransitionSet(epq.Element.Si,
                              epq.XRayTransitionSet.K_FAMILY)
siTrs = siTS.getTransitions()
for tr in siTrs:
    xrts.append(tr)


xp.update(mc3.configureXRayAccumulators(xrts,
                                        charAccum=True, 
                                        charFluorAccum=True, 
                                        bremFluorAccum=True))
# xp.update(mc3.configureEmissionImages(xrts,imgSzUm*1.0e-6, imgSzPx))
# xp.update(mc3.configurePhiRhoZ(imgSzUm*1.0e-6))
# xp.update(mc3.configureTrajectoryImage(imgSzUm*1.0e-6, imgSzPx))
# xp.update(mc3.configureVRML(nElectrons = 40))
xp.update(mc3.configureOutput(outDir))
xp.update(mc3.configureGun(beamSzNm*1.0e-9))

ts = time.time()

# first sim bare Si
spc = mc3.simulate(si, det, e0, dose, True, nTraj, True, True, xp)
spc = epq.SpectrumUtils.addNoiseToSpectrum(spc, 1.0)
spc = wrap(spc)

sName = "Si-%g-kV-%d-traj" % ( e0, nTraj)
spc.rename(sName)
res = anaConSi(spc, det, digits=2, display=False)
cI  = res["C"]
siI = res["Si"]
    
lT.append(0.0)
lPICmu.append(cI[0])
lPICuc.append(cI[1])
lPISimu.append(siI[0])
lPISiuc.append(siI[1])

te = time.time()

d = (te-ts) / 60.0

print("t = 0 nm, %.1f min" % (d) )

path = outDir + "/*/"

# clean up any un-needed MC* output directories
dirs = glob.glob(path)
for dir in dirs:
    dir = dir.replace('\\', '/')
    print(dir)
    shutil.rmtree(dir)





iCount += 1

lThickC = range( tNmStepC, tNmMaxC+tNmStepC, tNmStepC)

for tNmC in lThickC:
    ts = time.time()
    spc = mc3.coatedSubstrate(c, tNmC*1.0e-3, si, det, e0, True,
                              nTraj, dose, True, True, xp)
    spc = epq.SpectrumUtils.addNoiseToSpectrum(spc, 1.0)
    spc = wrap(spc)
    sName = "%g-nm-C-on-Si-%g-kV-%d-traj" % (tNmC, e0, nTraj)
    spc.rename(sName)
    res = anaConSi(spc, det, digits=2, display=False)
    cI  = res["C"]
    siI = res["Si"]
    lT.append(tNmC)
    lPICmu.append(cI[0])
    lPICuc.append(cI[1])
    lPISimu.append(siI[0])
    lPISiuc.append(siI[1])
    iCount += 1
    te = time.time()
    d = (te-ts) / 60.0
    print( "t = %g nm, %.1f min " % (tNmC, d) ) 
    # clean up any un-needed MC* output directories
    path = outDir + "/*/"
    dirs = glob.glob(path)
    for dir in dirs:
        dir = dir.replace('\\', '/')
        print(dir)
        shutil.rmtree(dir)



basFile ="C-ctd-Si-%g-kV-%g-Traj.csv" % (e0, nTraj)
strOutFile = outDir + "/" + basFile

f=open(strOutFile, 'w')
strLine = "t.C.nm, "
strLine = strLine +  "C.Int.mu,  C.Int.unc, "
strLine = strLine +  "Si.Int.mu, Si.Int.unc\n"

f.write(strLine)
for i in range(iCount):
    strLine = "%.3f" % lT[i] + ","
    strLine = strLine + "%.3f" % lPICmu[i] + ","
    strLine = strLine + "%.3f" % lPICuc[i] + ","
    strLine = strLine + "%.3f" % lPISimu[i] + ","
    strLine = strLine + "%.3f" % lPISiuc[i] + "\n"
    f.write(strLine)  
f.close()



# clean up cruft
shutil.rmtree(rptDir)
print "Done!"

end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg