# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

sim-K309-quant.py

A script to quantify K309 glass based on a video by N. Ritchie

  Date      Who   Ver    Comment
----------  ---  -----   ------------------------------------------
2019-12-16  JRM  0.0.1   Initial test. Required 38.168 min
                         Elapse: 0:38:10.4 min for
                         10,000 trajectories. This version requires
                         manual saving spectra. Added a small amount
                         of Ti to match the video...

1> tabulate(selected())
Name  O        Al      Si       Ca       Ti     
K309  37.9945  7.8947  18.8471  10.7828  0.0527
      Fe       Zn      Ba       Total
      10.5587  0.8665  13.9387  100.9356

For some reason the residual looks like Bremsstrahlung...

"""

__revision__ = "$Id: sim-K309-quant.py John R. Minter $"
__version__ = "0.0.1"

import sys
import os
import glob
import shutil
import fnmatch
import time
import math
import csv
import codecs

sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3, FluorescenceXRayGeneration3, XRayTransport3", None)
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.microanalysis.NISTMonte.Gen3 as nm3
import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQLibrary.Detector as epd
import gov.nist.microanalysis.Utility as epu
import gov.nist.microanalysis.EPQTools as ept
import java.io as jio
import java.lang as jl
import java.util as ju

from java.lang import Double


import dtsa2 as dt2
import gov.nist.microanalysis.dtsa2 as gdtsa2
import dtsa2.mcSimulate3 as mc3
import dtsa2.jmGen as jmg
import dtsa2.jmMC3 as jm3

det = findDetector("Oxford p4 05eV 4K")
e0       =    15    # kV
nTraj    = 10000    # trajectories
lt       =   500    # sec
pc       =     5.0 # nA
dose     = pc * lt  # na-sec"
bSaveSpc = True


gitHomeDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/sim-K309-quant/"
datDir = gitHomeDir + relPrj
print(gitHomeDir)
print(datDir)

"""
relPrj = "/Documents/work/dtsa2/sim-C-on-Si-7kV"
datDir = homDir + relPrj + "/msa-%g" % (nTraj)
jmg.ensureDir(datDir)
rptDir = homDir + relPrj + "/sim-C-on-Si-7kV Results"
csvFil = homDir + relPrj + "/dtsa2-C-on-Si-%g-kV-kratios-%g-traj.csv" % (e0, nTraj)
"""

start = time.time()

DataManager.clearSpectrumList()

"""
From N. Ritchie 2019-12-16

K309
Density 3 g/cm³

Element Z       Mass Frac       Norm Mass       Atom Frac
Oxygen          8       0.3872          0.3872          0.615282
Aluminum        13      0.0794          0.0794          0.0748163
Silicon         14      0.187           0.187           0.169278
Calcium         20      0.1072          0.1072          0.0680035
Iron            26      0.1049          0.1049          0.0477566
Barium          56      0.1343          0.1343          0.0248635
--      19.14   1       1.0     1.0

"""


l_comps = [epq.Element.O,  epq.Element.Al, epq.Element.Si,
           epq.Element.Ca, epq.Element.Fe, epq.Element.Ti,
           epq.Element.Ba]

# Note add a bit of Ti (impurity) in video not in list

l_mfs   = [0.3872,         0.0794,         0.187,
           0.1072,         0.1049,         0.0006,
           0.1343]

k309 = jmg.create_material_from_mf(l_comps, l_mfs, 2.6, "K309")

si = material("Si",density=2.32)
caf2 = material("CaF2", density=3.18)
fe = material("Fe", density=7.874)
ti = material("Ti", density=4.506)
zn = material("Zn", density=7.14)
baf2 = material("BaF2", density=4.89)
al2o3 = material("Al2O3", density=3.95)

xrts = []

trs = mc3.suggestTransitions(k309, e0)
for tr in trs:
    xrts.append(tr)

k309_spc = jm3.simBulkStd(k309, det, e0, nTraj, lt, pc, True)
k309_spc.display()

si_spc = jm3.simBulkStd(si, det, e0, nTraj, lt, pc, True)
si_spc.display()

caf2_spc = jm3.simBulkStd(caf2, det, e0, nTraj, lt, pc, True)
caf2_spc.display()

fe_spc = jm3.simBulkStd(fe, det, e0, nTraj, lt, pc, True)
fe_spc.display()

ti_spc = jm3.simBulkStd(ti, det, e0, nTraj, lt, pc, True)
ti_spc.display()

zn_spc = jm3.simBulkStd(zn, det, e0, nTraj, lt, pc, True)
zn_spc.display()

baf2_spc = jm3.simBulkStd(baf2, det, e0, nTraj, lt, pc, True)
baf2_spc.display()

al2o3_spc = jm3.simBulkStd(al2o3, det, e0, nTraj, lt, pc, True)
al2o3_spc.display()






print "Done!"

end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg




