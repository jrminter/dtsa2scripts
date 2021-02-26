# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

nist-k1526.py

A reproducible script to NIST K1526 glass spectrum

  Date      Who   Ver    Comment
----------  ---  -----   ------------------------------------------
2021-02-08  JRM  0.0.1   Initial test
                         For 10,000 traj
                         This script required 56.903 min
                         Elapse: 0:56:54.1
                         Simulate NIST K1526 glass spectrum

"""

__revision__ = "$Id: sim_310_steel_and_quant_spectra.py John R. Minter $"
__version__ = "0.0.1"

import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3,BremsstrahlungXRayGeneration3,FluorescenceXRayGeneration3, XRayTransport3", None)

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


start = time.time()

DataManager.clearSpectrumList()

def defineMat(elms,qty,name,density=None):
	c=epq.Composition(map(element,elms),qty,name)
	if density:
		c=epq.Material(c,epq.ToSI.gPerCC(density))
	return(c)


homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/sim-nist-k1526/"
e0      = 20.0
c_thick = 0.020 # 20 nm
n_traj  = 10000
dose    = 120.0
det     = findDetector("Oxford p4 05eV 4K")
ctg     = material("C", 1.8)


k1526 = defineMat(("Li",   "F",     "Na",    "Al",    "P",      "Br",    "O"      ),
	              (0.03908, 0.20274, 0.05475, 0.06135, 0.211245, 0.06441, 0.237405),
	              "K-15275", 3.30)

# compute spectrum
spc_k1526 = mc3.coatedSubstrate(ctg, c_thick, k1526, det, e0, True, n_traj, dose, True, True, {})

spc_k1526.rename("K-1526")
spc_k1526.setAsStandard(k1526)
spc_k1526.display()

sName = "20-nm-C-on-K1526-glass-at-%g-kV" % (e0)
spc_k1526.rename(sName)
spc_k1526.display()

end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg
