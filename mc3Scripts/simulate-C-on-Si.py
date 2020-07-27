# simulate-C-on-Si.py
#
#   Date      Who  Comment
# ----------  ---  -----------------------------------------------
# 2016-09-02  JRM  Simulate a thickness series of C on Si
#                  0 to 100 nm C in 5 nm steps
# This script required 17.066 min on ROCPW7ZC5C42 for 10K traj
# Elapse: 0:17:04.4

import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3,FluorescenceXRayGeneration3, XRayTransport3", None)


import os

import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQLibrary.Detector as epd
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.microanalysis.NISTMonte.Gen3 as nm3
import gov.nist.microanalysis.EPQTools as et
import gov.nist.microanalysis.dtsa2 as dt2

import dtsa2.jmGen as jmg
import dtsa2.mcSimulate3 as mc3
import dtsa2.hyperTools as ht

import shutil
import time

import java.util as jutil
import java.io as jio
import java.nio.charset as cs
import string

def anaCSi(spc, det, digits=2, display=True):
    """anaCSi

    Strip continuum and analyze a spectrum for C and Si

    Parameters
    ----------
    spc: A DTSA-II scriptable spectrum
        The spectrum to process
    det: The detector

    Returns
    -------
    A dictionary of tuples for C and Si where each tuple is the
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
    cInt = jmg.compPeakIntegral( bks, 0.282, 180, 1)
    siInt = jmg.compPeakIntegral(bks, 1.746, 240, 1)

    out = { "C": (round(cInt[0]/dose, digits), round(cInt[1]/dose, digits)),
            "Si": (round(siInt[0]/dose, digits), round(siInt[1]/dose, digits))
          }

    return out


start = time.time()

gitHom = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/mc3Scripts"
prjDir = gitHom + relPrj
datDir = prjDir + "/dat"
csvDir = datDir + "/csv"
jmg.ensureDir(datDir)
jmg.ensureDir(csvDir)

os.chdir(prjDir)
rptDir = prjDir + '/simulate-C-on-Si Results/'


det      = findDetector("Oxford p4 05eV 2K")
e0       =     5    # kV
nTraj    = 10000    # trajectories
lt       =   100    # sec
pc       =     0.25 # nA
tNmStepC =    10
tNmCMax  =   100



dose = pc * lt  # na-sec"

DataManager.clearSpectrumList()

c = material("C", density=2.2)
si = material("Si", density=2.3296)

lThickC = range(0, tNmCMax+tNmStepC, tNmStepC)

lPICmu = [] # an array for mean C peak integral
lPICuc = [] # an array for C peak integral uncertainty
lPISimu = [] # an array for mean Si peak integral
lPISiuc = [] # an array for Si peak integral uncertainty


iCount = 0

for tNmC in lThickC:
    if tNmC == 0:
        layers = [ [si, 50.0e-6] ]
    
    else:
        layers = [ [c, tNmC*1.0e-9],
                   [si, 50.0e-6]
                 ]

    spc = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                        dose=dose, sf=True, bf=True, xtraParams={})
    
    if tNmC < 1:
        sName = "Si-%g-kV" % (e0)
        spcNa = "Si"
    else:
        sName = "%g-nm-C-on-Si-%g-kV" % (tNmC, e0)
        spcNa = "%g nm C on Si" % (tNmC)

    spc.rename(spcNa)
    res = anaCSi(spc, det, digits=2, display=True)
    cI = res["C"]
    siI = res["Si"]
    lPICmu.append(cI[0])
    lPICuc.append(cI[1])
    lPISimu.append(siI[0])
    lPISiuc.append(siI[1])
    iCount += 1
    print(res)

basFile ="C-on-Si-%gkV-%g-Traj.csv" % (e0, nTraj)
strOutFile = csvDir + "/" + basFile

f=open(strOutFile, 'w')
strLine = "t.C.nm, "
strLine = strLine +  "C.Int.mu, C.Int.unc, "
strLine = strLine +  "Si.Int.mu, Si.Int.unc\n"
f.write(strLine)
for i in range(iCount):
    strLine = "%.3f" % lThickC[i] + ","
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
