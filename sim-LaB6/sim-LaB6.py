# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

sim-LaB6.py

A reproducible script to simulate a LaB6 and LaPO4 spectra

  Date      Who   Ver    Comment
----------  ---  -----   ------------------------------------------
2021-04-29  JRM  0.0.1   Initial test
                         For 10,000 traj
                         This script required min
                         Elapse: 0:56:54.1


"""

__revision__ = "$Id: sim-LaB6.py John R. Minter $"
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
wrkDir = homDir + "/Documents/git/dtsa2Scripts/sim-LaB6/"
e0      = 20.0
c_thick = 0.020 # 20 nm
n_traj  = 10000
dose    = 120.0
det     = findDetector("Oxford p4 05eV 4K")
ctg     = material("C", 1.8)


lab6  = defineMat(("La"   ,  "B",    ),
	              (0.68168, 0.31832, ),
	               "LaB6", 4.72)

lapo4 = defineMat(("La"   ,     "P",  "O"    ),
	              (0.59393, 0.13243,  0.27364,),
	              "LaPO4", 4.93)
# compute spectrum
spc_lab6 = mc3.coatedSubstrate(ctg, c_thick, lab6, det, e0, True, n_traj, dose, True, True, {})
spc_lab6.rename("LaB6")
# spc_lab6.setAsStandard(lab6)
spc_lab6.display()

sName = "20-nm-C-on-LaB6-%g-kV" % (e0)
spc_lab6.rename(sName)
spc_lab6.display()

spc_lapo4 = mc3.coatedSubstrate(ctg, c_thick, lapo4, det, e0, True, n_traj, dose, True, True, {})
spc_lapo4.rename("LaPO4")
# spc_lapo4.setAsStandard(lapo4)
spc_lapo4.display()

sName = "20-nm-C-on-LaB6-%g-kV" % (e0)
spc_lapo4.rename(sName)
spc_lapo4.display()

end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg
