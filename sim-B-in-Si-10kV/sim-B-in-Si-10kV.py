# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

sim-B-in-Si-10kV.py

  Date      Who   Ver    Comment
----------  ---  -----   ------------------------------------------
2020-07-05  JRM  0.0.1   Initial test

"""

__revision__ = "$Id: sim-B-in-Si-10kV.py John R. Minter $"
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

e0      =    10.0
n_traj  = 10000
dose    =   120.0
det     =  findDetector("Oxford p4 05eV 2K")

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/sim-B-in-Si-10kV/"

def defineMat(elms,qty,name,density=None):
  c=epq.Composition(map(element,elms),qty,name)
  if density:
    c=epq.Material(c,epq.ToSI.gPerCC(density))
  return(c)

DataManager.clearSpectrumList()

"""
Composition
"""
bsi   = defineMat(("B", "Si", ), (0.500000, 0.500000), "BSi",   2.66)
b2si  = defineMat(("B", "Si", ), (0.565020, 0.434980), "B2Si",  2.66)
b5si  = defineMat(("B", "Si", ), (0.341924, 0.658076), "B5Si",  2.66)
b7si  = defineMat(("B", "Si", ), (0.270675, 0.729326), "B7Si",  2.66)
b95si = defineMat(("B", "Si", ), (0.026618, 0.973382), "B95Si", 2.66)



bsi_spc  = simulate(bsi, det, e0, 60.0, withPoisson=False)
bsi_spc.rename("BSi  ")
bsi_spc.display()

b2si_spc  = simulate(b2si, det, e0, 60.0, withPoisson=False)
b2si_spc.rename("B2Si ")
b2si_spc.display()

b5si_spc  = simulate(b5si, det, e0, 60.0, withPoisson=False)
b5si_spc.rename("B5Si ")
b5si_spc.display()

b7si_spc  = simulate(b7si, det, e0, 60.0, withPoisson=False)
b7si_spc.rename("B7Si ")
b7si_spc.display()


b95si_spc = simulate(b95si, det, e0, 60.0, withPoisson=False)
b95si_spc.rename("B95Si")
b95si_spc.display()





end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg





