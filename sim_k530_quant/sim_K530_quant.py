# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

sim_K530_and_quant.py

A reproducible script to compute spectra and quantify K530 glass
and compare to K412

  Date      Who   Ver    Comment
----------  ---  -----   ------------------------------------------
2020-06-17  JRM  0.0.1   Initial test For 10,000 traj

"""

__revision__ = "$Id: sim_K530_glass_spectra.py John R. Minter $"
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

e0      = 20.0
c_thick = 0.020 # 20 nm
n_traj  = 10000
dose    = 120.0
det     = findDetector("Oxford p4 05eV 4K")

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/sim_K530_quant"

ctg     = material("C", 1.8)

"""
Composition of K-412 from Nicholas Ritchie
"""
k412nr  = defineMat(("O",    "Mg",   "Al",  "Si",   "Ca",   "Fe", ),
	                (0.4312, 0.1175, 0.0495, 0.2138, 0.1099, 0.0781),
                    "K412-NR", 2.66)

"""
Composition of K-412 from John Fournelle 
"""


k412jf = defineMat(("O",   "Mg",   "Al",  "Si",  "Ca",  "Fe", ),
                   (0.4360,0.1166,0.0491,0.2151,0.1090,0.0774),
	                "K412-JF", 2.66)
"""
Composition of K-530 from John Fournelle 
"""

k530jf  = defineMat(("O",   "Mg",   "Al",  "Si",  "Ca",  "Fe", ),
	                 (0.4318,0.1152,0.0510,0.2151,0.1068,0.0801),
	                "K530-JF", 2.66)
"""
The CalcZAF Standard program doe the computations for us in mf and oxide

Average Total Oxygen:         .432     Average Total Weight%:    1.000
Average Calculated Oxygen:    .432     Average Atomic Number:   12.730
Average Excess Oxygen:        .000     Average Atomic Weight:   22.039

ELEM:        O     MgO   Al2O3    SiO2     CaO     FeO
XRAY:      ka      ka      ka      ka      ka      ka 
OXWT:     .000    .191    .096    .460    .149    .103
ELWT:     .432    .115    .051    .215    .107    .080
KFAC:    .0017   .0008   .0003   .0016   .0010   .0007
ZCOR:   2.5028  1.4996  1.4628  1.3079  1.0816  1.1834
AT% :   59.477  10.446   4.166  16.879   5.873   3.161
24 O:   24.000   4.215   1.681   6.811   2.370   1.276



"""

# simulate K412-NR
spc_k412nr = mc3.coatedSubstrate(ctg, c_thick, k412nr, det, e0,
                                 True, n_traj, dose, True, True, {})
sName = "%g nm C on K412-NR at %g kV" % (1000*c_thick, e0)
spc_k412nr.rename(sName)
spc_k412nr.display()
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_k412nr.setAsStandard(k412nr)
spc_k412nr.save(fi)


# simulate K412-JF
spc_k412jf = mc3.coatedSubstrate(ctg, c_thick, k412jf, det, e0,
                                 True, n_traj, dose, True, True, {})
sName = "%g nm C on K412-JF at %g kV" % (1000*c_thick, e0)
spc_k412jf.rename(sName)
spc_k412jf.setAsStandard(k412jf)
spc_k412jf.display()
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_k412jf.save(fi)

# simulate K530-JF
spc_k530jf = mc3.coatedSubstrate(ctg, c_thick, k530jf, det, e0,
                                 True, n_traj, dose, True, True, {})
sName = "%g nm C on K530-JF at %g kV" % (1000*c_thick, e0)
spc_k530jf.rename(sName)
spc_k530jf.setAsStandard(k530jf)
spc_k530jf.display()
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_k530jf.save(fi)


end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg





