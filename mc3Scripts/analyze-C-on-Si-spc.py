# analyze-C-on-Si-spc.py
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
import glob

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
relPrj = "/dtsa2Scripts"
prjDir = gitHom + relPrj
wrkDir = prjDir + "/mc3Scripts"
spcDir = prjDir + "/spc/Sirion/Oxford-P4-05eV-2K/Si-Cal-5kV"
csvDir = wrkDir + "/dat/csv"
jmg.ensureDir(wrkDir)

os.chdir(wrkDir)

rptDir = wrkDir + '/analyze-C-on-Si-spc Results/'


det = findDetector("Oxford p4 05eV 2K")
e0  = 5    # kV

DataManager.clearSpectrumList()

lPICmu = [] # an array for mean C peak integral
lPICuc = [] # an array for C peak integral uncertainty
lPISimu = [] # an array for mean Si peak integral
lPISiuc = [] # an array for Si peak integral uncertainty

print(spcDir + '/*.msa')
for name in glob.glob(spcDir + '/*.msa'):
    name = name.replace('\\', '/')
    bn = os.path.basename(name)
    na = bn.split('.msa')[0]
    print(na)
    spc = wrap(readSpectrum(name))
    sp = spc.getProperties()
    pc = sp.getNumericProperty(epq.SpectrumProperties.FaradayBegin)
    lt = sp.getNumericProperty(epq.SpectrumProperties.LiveTime)
    spc.display()
    res = anaCSi(spc, det, digits=2, display=True)
    cI = res["C"]
    siI = res["Si"]
    ratio = cI[0]/siI[0]
    print((na, pc, siI[0], siI[1], ratio))

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
