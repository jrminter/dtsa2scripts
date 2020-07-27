# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

sim_azo_and_quant_spectra.py

A reproducible script to compute spectra and quantify AZO2 and AZO5

  Date      Who   Ver    Comment
----------  ---  -----   ------------------------------------------
2020-06-10  JRM  0.0.1   Initial test
                         For 10,000 traj
                         This script required 35.143 min
                         Elapse: 0:35:08.6

"""

__revision__ = "$Id: sim_k411_and_k412_quant_spectra.py John R. Minter $"
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
wrkDir = homDir + "/Documents/git/dtsa2Scripts/sim_azo_and_quant_spectra/"



ctg = material("C", 1.8)

al = material("Al", 2.7)

azo2 = defineMat(("O","Al","Zn",),
	               (0.1984,0.0067,0.7949,),
	               "AZO2", 5.61)
azo5 = defineMat(("O","Al","Zn",),
                 (0.2013,0.0170,0.7917,),
                 "AZO5", 5.61)
zno = defineMat(("O","Zn",),
                 (0.1965,0.8035,),
                 "ZnO", 5.61)

"""
coatedSubstrate(coating, thickness, substrate, det, e0=20.0,
                withPoisson=True, nTraj=defaultNumTraj,
                dose=defaultDose, sf=defaultCharFluor,
                bf=defaultBremFluor, xtraParams=defaultXtraParams)
"""

spc_al = mc3.coatedSubstrate(ctg, c_thick, al, det, e0,
                               True, n_traj, dose, True, True, {})
sName = "%g nm C on Al at %g kV" % (1000*c_thick, e0)
spc_al.rename(sName)
spc_al.display()
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_al.save(fi)

spc_azo2 = mc3.coatedSubstrate(ctg, c_thick, azo2, det, e0,
                               True, n_traj, dose, True, True, {})
sName = "%g nm C on AZO2 at %g kV" % (1000*c_thick, e0)
spc_azo2.rename(sName)
spc_azo2.display()
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_azo2.save(fi)

spc_azo5 = mc3.coatedSubstrate(ctg, c_thick, azo5, det, e0,
                               True, n_traj, dose, True, True, {})
sName = "%g nm C on AZO5 at %g kV" % (1000*c_thick, e0)
spc_azo5.rename(sName)
spc_azo5.display()
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_azo5.save(fi)

spc_zno = mc3.coatedSubstrate(ctg, c_thick, zno, det, e0,
                               True, n_traj, dose, True, True, {})
sName = "%g nm C on ZnO at %g kV" % (1000*c_thick, e0)
spc_zno.rename(sName)
spc_zno.display()
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_zno.setAsStandard(zno)
spc_zno.save(fi)


end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg

