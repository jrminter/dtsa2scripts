# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

sim-pyrex.py

A reproducible script to simulate a Pyrex spectrum

  Date      Who   Ver    Comment
----------  ---  -----   ------------------------------------------
2021-06-09  JRM  0.0.1   Initial test
                         For 10000 traj on MacOS orig laptop
                         This script required 17.327 min
                         Elapse: 0:17:19.7
                         Running /Users/jrminter/Documents/git/dtsa2Scripts/sim-Pyrex/sim-pyrex.py
"""

__revision__ = "$Id: sim-pyrex.py John R. Minter $"
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
wrkDir = homDir + "/Documents/git/dtsa2Scripts/sim-pyrex/"
e0      = 20.0
c_thick = 0.020 # 20 nm
n_traj  = 10000
dose    = 120.0
det     = findDetector("Oxford p4 05eV 4K")
ctg     = material("C", 1.8)


pyrex   = defineMat((   "B", "O",      "Na",   "Mg",   "Al",   "Si",
                        "Cl", "Ca", "Fe"   ),
                    (0.0392, 0.5384, 0.0312, 0.0003, 0.0117, 0.3772,
                     0.0010, 0.0007, 0.0003),
                     "pyrex", 2.5)

# compute spectrum
spc_pyrex = mc3.coatedSubstrate(ctg, c_thick, pyrex, det, e0, True,
                                n_traj, dose, True, True, {})

spc_pyrex.rename("pyrex")
spc_pyrex.setAsStandard(pyrex)
spc_pyrex.display()

sName = "20-nm-C-on-Pyrex-at-%g-kV" % (e0)
spc_pyrex.rename(sName)
spc_pyrex.display()

end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg
