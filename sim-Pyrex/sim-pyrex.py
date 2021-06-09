# -*- coding: utf-8 -*-
<<<<<<< HEAD
#         1         2         3         4         5         6         7 |
# 23456789012345678901234567890123456789012345678901234567890123456789012
#
# sim-pyrex.py
# jrm 2018-09-13  starting a run with 50,000 traj
#                This script required 379.140 min
#                Elapse: 6:19:09.2
import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3, FluorescenceXRayGeneration3, XRayTransport3", None)
=======

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
import dtsa2.mcSimulate3 as mc3
import dtsa2.jmGen as jmg
import dtsa2.jmMC3 as jm3
import java.util as jutil
import java.io as jio
import java.nio.charset as cs
import os
import shutil
import time
import datetime


det = findDetector("Oxford p4 05eV 2K")
e0       =     5    # kV
nTraj    =  1000    # trajectories
lt       =   500    # sec
pc       =     5.0 # nA
dose     = pc * lt  # na-sec"
bSaveSpc = True

homDir = os.environ['HOME']
relPrj = "/Documents/work/dtsa2/sim-C-on-Pyrex-5kV"
datDir = homDir + relPrj + "/msa-%g" % (nTraj)
jmg.ensureDir(datDir)

def sim_amc_coated_mat(mat, det, e0, nTraj, lt=100, pc=1.0, tc=20.0):
    """sim_amc_coated_mat(mat, det, e0, nTraj, lt=100, pc=1.0, tc=20.0)

    Use mc3 multilayer simulation to simulate an am-C-ctd specimen

    Parameters
    ----------
    mat - a dtsa material.
        Note the material must have an associated density. It should have a useful name.
    det - a dtsa detector
        Here is where we get the detector properties and calibration
    e0 - float
        The accelerating voltage in kV
    nTraj - integer
        The number of trajectories to run
    lt - integer (100)
        The live time (sec)
    pc - float (1.0)
        The probe current in nA
    tc - float (20.0)
        C thickness in nm


    Returns
    -------
    sim - DTSA scriptable spectrum 
        The simulated spectrum
    
    Example
    -------
    import dtsa2.jmMC3 as jm3
    det = findDetector("Oxford p4 05eV 2K")
    si = material("Si", density=2.3296)
    a = jm3.simCarbonCoatedStd(mgo, det, 20.0, 100, 100, 1.0, 20.0)
    a.display()

    """
    dose = pc * lt  # na-sec"
    amc = material("C", density=1.35)
    
    amcThickComment = "amC Thickness = %g nm %g trajectories" % (tc, nTraj)
    layers = [ [amc, tc*1.0e-9],
               [mat, 50.0e-6]
             ]
    xrts = []

    trs = mc3.suggestTransitions(amc, e0)
    for tr in trs:
       xrts.append(tr)

    trs = mc3.suggestTransitions(mat, e0)
    for tr in trs:
        xrts.append(tr)

    xtraParams={}
    xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))

    sim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                        dose=dose, sf=True, bf=True, xtraParams=xtraParams)
    sName = "%g-nm-amC-on-%s-%g-kV" % (tc, mat, e0)
    sim.rename(sName)
    sim.getProperties().setTextProperty(epq.SpectrumProperties.SpectrumComment, amcThickComment)
    return sim
=======
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

>>>>>>> 157806ec9007eba590a295c404035244eb5759bc

start = time.time()

DataManager.clearSpectrumList()

<<<<<<< HEAD
c   = material("C", density=2.2)
si  = material("Si", density=2.329)



=======
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
>>>>>>> 157806ec9007eba590a295c404035244eb5759bc
