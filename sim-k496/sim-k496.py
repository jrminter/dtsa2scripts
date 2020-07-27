# test-sim.py
#   Date      Who  Comment
# ----------  ---  -----------------------------------------------
# 2020-06-25  JRM  Test simulation of K497 quant spectra
#                  Modified to store the spectra in "spc" dir
#
# 2020-06-28  JRM  Ran again after recomputing mass fractions for
#                  K496 and K497. Also added these to the DTSA2 database
# 
# This script required 76.014 min for 10,000 trajectories
# on jrmSimulation PC 16GB RAM  HP Obelisk 
# This script required 74.196 min
# ...or 1.237 hr
# Elapse: 1:14:11.7
# Spectra are in on My C drive in \Documents\git\dtsa2Scripts\mc3Scripts\dat\spc

import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3,FluorescenceXRayGeneration3, XRayTransport3", None)


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

gitHom = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts"
prjDir = gitHom + relPrj
wrkDir = prjDir + "/sim-k496"
jmg.ensureDir(wrkDir)
os.chdir(wrkDir)
# rptDir = wrkDir + '/Test-sim Results/'
# Save spectra in spcDir directory
spcDir = wrkDir + "/spc"
jmg.ensureDir(spcDir)

det     = findDetector("Oxford p4 05eV 4K")
e0      = 15      # kV
nTraj   = 100  # trajectories
lt      = 60      # sec
pc      = 0.25    # nA
tNmC    = 20



dose = pc * lt  # na-sec"

DataManager.clearSpectrumList()

ctg     = material("C", 1.8)

# K 496

"""
Al  0.06470
Mg  0.06650
P   0.32980
O   0.05390
==  ========
Tot 1.0

N.B. John Donovan notes:
From NBS (NIST), Microanalysis glasses
2 x 2 x 20 mm rod
Traces by EPMA (NBS) by R. Marinenko
**Oxygen by difference**


"""

l_comps = [epq.Element.Al,  epq.Element.Mg, epq.Element.P, epq.Element.O]
l_mfs   = [ 0.06470,        0.06650,        0.32980,       0.5390]
k496 = jmg.create_material_from_mf(l_comps, l_mfs,  2.6, "K-496")



# K497

l_comps = [epq.Element.Al,  epq.Element.Mg, epq.Element.P,  epq.Element.O,
           epq.Element.Pb,  epq.Element.Si, epq.Element.Ti, epq.Element.Fe,
           epq.Element.Zr,  epq.Element.Ce, epq.Element.Ta]
l_mfs   = [ 0.045950,       0.055450,       0.211790,        0.680870,
            0.000860,       0.000960,       0.000950,        0.000970,
            0.000730,       0.001390,       0.000080]

"""
Al  0.045950
Mg  0.055450
P   0.211790
O   0.680870
Pb  0.000860
Si  0.000960
Ti  0.000950
Fe  0.000970
Zr  0.000730
Ce  0.001390
Ta  0.000080
---  --------
Tot  1.000000
"""


k497 = jmg.create_material_from_mf(l_comps, l_mfs,  2.6, "K-497")

# sio2 = material("SiO2", density=2.65)
c = material("C", density=2.1)
# aupd = epq.Material(epq.Composition([epq.Element.Au, epq.Element.Pd],
#                                    [0.60,0.40]),
#                                    epq.ToSI.gPerCC(0.6*19.30+0.4*11.9))
#

# K496 Sim
#
layers = [ [c, 2.0e-9],
           [k496, 50.0e-6]
         ]

# multiFilm(layers, det, e0=20.0, withPoisson=True,
#           nTraj=defaultNumTraj, dose=defaultDose,
#           sf=defaultCharFluor, bf=defaultBremFluor,
#            xtraParams=defaultXtraParams)



k496Sim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                        dose=dose, sf=True, bf=True, xtraParams={})
sName = "%g-nm-C-on-K496" % tNmC
k496Sim.rename(sName)
k496Sim.setAsStandard(k496)
k496Sim.display()
fi =  spcDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (nTraj)
k496Sim.save(fi)


# K-497 Sim
layers = [ [c, 2.0e-9],
           [k497, 50.0e-6]
         ]

k497Sim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                        dose=dose, sf=True, bf=True, xtraParams={})
sName = "%g-nm-C-on-K497" % tNmC
k497Sim.rename(sName)
k497Sim.setAsStandard(k497)
k497Sim.display()
fi =  spcDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (nTraj)
k497Sim.save(fi)

# Al2O3 Sim
al2o3 = material("Al2O3", 3.95)
layers = [ [c, 2.0e-9],
           [al2o3, 50.0e-6]
         ]
xrts = mc3.suggestTransitions(al2o3)
al2o3Sim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                         dose=dose, sf=True, bf=True, xtraParams={})
sName = "%g-nm-C-on-Al2O3" % tNmC
al2o3Sim.rename(sName)
al2o3Sim.setAsStandard(al2o3)
al2o3Sim.display()
fi =  spcDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (nTraj)
al2o3Sim.save(fi)

# MgO Sim
mgo   = material("MgO", 3.58)
xrts = mc3.suggestTransitions(mgo)
layers = [ [c, 2.0e-9],
           [mgo, 50.0e-6]
         ]

mgoSim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                         dose=dose, sf=True, bf=True, xtraParams={})
sName = "%g-nm-C-on-MgO" % tNmC
mgoSim.rename(sName)
mgoSim.setAsStandard(mgo)
mgoSim.display()
fi =  spcDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (nTraj)
mgoSim.save(fi)

# PBS Sim
pbs = material("PbS", 7.6)
xrts = mc3.suggestTransitions(pbs)
layers = [ [c, 2.0e-9],
           [pbs, 50.0e-6]
         ]
pbsSim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                       dose=dose, sf=True, bf=True, xtraParams={})
sName = "%g-nm-C-on-PbS" % tNmC
pbsSim.rename(sName)
pbsSim.setAsStandard(pbs)
pbsSim.display()
fi =  spcDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (nTraj)
pbsSim.save(fi)

# Si Sim
si = material("Si", 2.32)
xrts = mc3.suggestTransitions(si)
layers = [ [c, 2.0e-9],
           [si, 50.0e-6]
         ]
siSim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                       dose=dose, sf=True, bf=True, xtraParams={})
sName = "%g-nm-C-on-Si" % tNmC
siSim.rename(sName)
siSim.setAsStandard(si)
siSim.display()
fi =  spcDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (nTraj)
siSim.save(fi)

# ZrO2 Sim
zro2 = material("ZrO2",  5.68)
xrts = mc3.suggestTransitions(zro2)
layers = [ [c, 2.0e-9],
           [zro2, 50.0e-6]
         ]
zro2Sim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                        dose=dose, sf=True, bf=True, xtraParams={})
sName = "%g-nm-C-on-ZrO2" % tNmC
zro2Sim.rename(sName)
zro2Sim.setAsStandard(zro2)
zro2Sim.display()
fi =  spcDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (nTraj)
zro2Sim.save(fi)

# TiO2 Sim
tio2 = material("TiO2", 4.23)
xrts = mc3.suggestTransitions(tio2)
layers = [ [c, 2.0e-9],
           [tio2, 50.0e-6]
         ]
tio2Sim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                        dose=dose, sf=True, bf=True, xtraParams={})
sName = "%g-nm-C-on-TiO2" % tNmC
tio2Sim.rename(sName)
tio2Sim.setAsStandard(tio2)
tio2Sim.display()
fi =  spcDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (nTraj)
tio2Sim.save(fi)


# Fe2O3 Sim
fe2o3 = material("Fe2O3", 5.24)
xrts = mc3.suggestTransitions(fe2o3)
layers = [ [c, 2.0e-9],
           [fe2o3, 50.0e-6]
         ]
fe2o3Sim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                         dose=dose, sf=True, bf=True, xtraParams={})
sName = "%g-nm-C-on-Fe2O3" % tNmC
fe2o3Sim.rename(sName)
fe2o3Sim.setAsStandard(fe2o3)
fe2o3Sim.display()
fi =  spcDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (nTraj)
fe2o3Sim.save(fi)

# Cs2O Sim
cs2o = material("Cs2O", 4.65)
xrts = mc3.suggestTransitions(cs2o)
layers = [ [c, 2.0e-9],
           [cs2o, 50.0e-6]
         ]
cs2oSim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                        dose=dose, sf=True, bf=True, xtraParams={})
sName = "%g-nm-C-on-Cs2O" % tNmC
cs2oSim.rename(sName)
cs2oSim.setAsStandard(cs2o)
cs2oSim.display()
fi =  spcDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (nTraj)
cs2oSim.save(fi)

# Ta Sim
ta = material("Ta", 16.65)
xrts = mc3.suggestTransitions(ta)
layers = [ [c, 2.0e-9],
           [ta, 50.0e-6]
         ]
taSim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                      dose=dose, sf=True, bf=True, xtraParams={})
sName = "%g-nm-C-on-Ta" % tNmC
taSim.rename(sName)
taSim.setAsStandard(ta)
taSim.display()
fi =  spcDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (nTraj)
taSim.save(fi)

print(fi)


# clean up cruft
# shutil.rmtree(rptDir)
print "Done!"

end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg

