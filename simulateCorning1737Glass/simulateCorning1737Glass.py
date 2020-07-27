# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

sim_k411_add_k412_quant_spectra.py

A reproducible script to compute spectrum for Corning 1737 Glass

  Date      Who   Ver    Comment
----------  ---  -----   ------------------------------------------
2020-06-17  JRM  0.0.1   Initial test
                         For 10,000 traj
"""

__revision__ = "$Id: simulateCorning1737Glass.py John R. Minter $"
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
  # Note: epq.Composition is designed to work with **mass fractions**
  # Note: epq.Composition(Element[] elms, double[] massFracs)
	c=epq.Composition(map(element,elms),qty,name)
	if density:
		c=epq.Material(c,epq.ToSI.gPerCC(density))
	return(c)

e0      = 20.0
c_thick = 0.020 # 20 nm
n_traj  = 100
dose    = 120.0
det     = findDetector("Oxford p4 05eV 4K")



homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/simulateCorning1737Glass/"



ctg = material("C", 1.8)

c1737 = defineMat(("Ba",  "Sr",  "As",  "Ca",   "Si",   "Al",   "Mg",  "O",    "B"),
	              (0.0863,0.1502,0.0428,0.02862,0.27676,0.08863,0.00486,0.47299,0.02254),
	               "Corning 1737 Glass",2.540)

"""
coatedSubstrate(coating, thickness, substrate, det, e0=20.0,
                withPoisson=True, nTraj=defaultNumTraj,
                dose=defaultDose, sf=defaultCharFluor,
                bf=defaultBremFluor, xtraParams=defaultXtraParams)
"""

spc_c1737 = mc3.coatedSubstrate(ctg, c_thick, c1737, det, e0,
                                True, n_traj, dose, True, True, {})
sName = "%g nm C on Corning 1737 Glass at %g kV" % (1000*c_thick, e0)
spc_c1737.rename(sName)
spc_c1737.display()
spc_c1737.setAsStandard(c1737)
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_c1737.save(fi)

end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg
