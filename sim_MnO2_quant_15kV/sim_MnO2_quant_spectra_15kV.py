# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

sim_MnO2_quant_spectra_15kV.py

A reproducible script to compute spectra and quantify MnO2 at 15 kV

  Date      Who   Ver    Comment
----------  ---  -----   ------------------------------------------
2020-06-11  JRM  0.0.1   Initial test
                         For 10,000 traj
                         This script required 21.286 min
                         Elapse: 0:21:17.1 MacOS

"""

__revision__ = "$Id: sim_MnO2_quant_spectra15kV.py John R. Minter $"
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
wrkDir = homDir + "/Documents/git/dtsa2Scripts/sim_MnO2_quant_15kV"


ctg = material("C", 1.8)

mno2 = material("MnO2", 5.03)

mn = material("Mn", 7.42)

sio2 = material("SiO2", 2.65)


"""
coatedSubstrate(coating, thickness, substrate, det, e0=20.0,
                withPoisson=True, nTraj=defaultNumTraj,
                dose=defaultDose, sf=defaultCharFluor,
                bf=defaultBremFluor, xtraParams=defaultXtraParams)
"""


spc_mno2 = mc3.coatedSubstrate(ctg, c_thick, mno2, det, e0,
                               True, n_traj, dose, True, True, {})
sName = "%g-nm-C-on-MnO2-at-%g-kV" % (1000*c_thick, e0)
spc_mno2.rename(sName)
spc_mno2.display()
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_mno2.setAsStandard(mno2)
print(fi)
spc_mno2.save(fi)

spc_sio2 = mc3.coatedSubstrate(ctg, c_thick, sio2, det, e0,
                               True, n_traj, dose, True, True, {})
sName = "%g-nm-C-on-SiO2-at-%g-kV" % (1000*c_thick, e0)
spc_sio2.rename(sName)
spc_sio2.display()
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_sio2.setAsStandard(sio2)
print(fi)
spc_sio2.save(fi)

spc_mn = mc3.coatedSubstrate(ctg, c_thick, mn, det, e0,
                               True, n_traj, dose, True, True, {})
sName = "%g-nm-C-on-Mn-at-%g-kV" % (1000*c_thick, e0)
spc_mn.rename(sName)
spc_mn.display()
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_mn.setAsStandard(mn)
print(fi)
spc_mn.save(fi)


end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg

