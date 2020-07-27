# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

sim_k411_add_k412_quant_spectra.py

A reproducible script to compute spectra and quantify K411 and K412 glasses

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
wrkDir = homDir + "/Documents/git/dtsa2Scripts/sim_k411_and_k412_quant_spectra/"



ctg = material("C", 1.8)

k411 = defineMat(("O","Mg","Si","Ca","Fe",),
	             (0.428600,0.089500,0.256700,0.111800,0.113400,),
	             "K411",2.600000)
k412 = defineMat(("O","Mg","Al","Si","Ca","Fe",),
	             (0.431202,0.117554,0.049477,0.213777,0.109914,0.078075,),
	             "K412",2.660000)
al2o3 = material("Al2O3", 3.95)

caf2 = material("CaF2", 3.18)

fe = material("Fe", 7.874)

si = material("Si", 2.329)

mgo = material("MgO", 3.58)

"""
coatedSubstrate(coating, thickness, substrate, det, e0=20.0,
                withPoisson=True, nTraj=defaultNumTraj,
                dose=defaultDose, sf=defaultCharFluor,
                bf=defaultBremFluor, xtraParams=defaultXtraParams)
"""

spc_k411 = mc3.coatedSubstrate(ctg, c_thick, k411, det, e0,
                               True, n_traj, dose, True, True, {})
sName = "%g nm C on K411 at %g kV" % (1000*c_thick, e0)
spc_k411.rename(sName)
spc_k411.display()
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_k411.save(fi)

spc_k412 = mc3.coatedSubstrate(ctg, c_thick, k412, det, e0,
                               True, n_traj, dose, True, True, {})
sName = "%g nm C on K412 at %g kV" % (1000*c_thick, e0)
spc_k412.rename(sName)
spc_k412.display()
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_k412.save(fi)

spc_al2o3 = mc3.coatedSubstrate(ctg, c_thick, al2o3, det, e0,
                               True, n_traj, dose, True, True, {})
sName = "%g nm C on Al2O3 at %g kV" % (1000*c_thick, e0)
spc_al2o3.rename(sName)
spc_al2o3.display()
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_al2o3.setAsStandard(al2o3)
spc_al2o3.save(fi)

spc_caf2 = mc3.coatedSubstrate(ctg, c_thick, caf2, det, e0,
                               True, n_traj, dose, True, True, {})
sName = "%g nm C on CaF2 at %g kV" % (1000*c_thick, e0)
spc_caf2.rename(sName)
spc_caf2.display()
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_caf2.setAsStandard(caf2)
spc_caf2.save(fi)

spc_fe = mc3.coatedSubstrate(ctg, c_thick, fe, det, e0,
                               True, n_traj, dose, True, True, {})
sName = "%g nm C on Fe at %g kV" % (1000*c_thick, e0)
spc_fe.rename(sName)
spc_fe.display()
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_fe.setAsStandard(fe)
spc_fe.save(fi)

spc_mgo = mc3.coatedSubstrate(ctg, c_thick, mgo, det, e0,
                               True, n_traj, dose, True, True, {})
sName = "%g nm C on MgO at %g kV" % (1000*c_thick, e0)
spc_mgo.rename(sName)
spc_mgo.display()
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_mgo.setAsStandard(mgo)
spc_mgo.save(fi)

spc_si = mc3.coatedSubstrate(ctg, c_thick, si, det, e0,
                               True, n_traj, dose, True, True, {})
sName = "%g nm C on Si at %g kV" % (1000*c_thick, e0)
spc_si.rename(sName)
spc_si.display()
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_si.setAsStandard(si)
spc_si.save(fi)




end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg

