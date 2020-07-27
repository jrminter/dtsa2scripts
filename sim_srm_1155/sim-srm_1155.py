# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

sim-srm_1155.py

A script to quantify K309 glass based on a video by N. Ritchie

  Date      Who   Ver    Comment
----------  ---  -----   ------------------------------------------
2019-12-25  JRM  0.0.1   Initial test.

ELEM:        C      Mn       P       S      Si      Cu      Ni      Cr
XRAY:      ka      ka      ka      ka      ka      ka      ka      ka
ELWT:     .046   1.630    .020    .018    .502    .169  12.180  18.450
KFAC:    .0001   .0163   .0002   .0002   .0033   .0016   .1166   .2128
ZCOR:   4.3792   .9981  1.3232  1.1449  1.5146  1.0843  1.0446   .8669
AT% :     .213   1.649    .036    .031    .993    .148  11.531  19.723

ELEM:        V      Mo      Co      Pb      Fe
XRAY:      ka      la      ka      ma      ka
ELWT:     .047   2.380    .101    .001  64.456
KFAC:    .0005   .0198   .0010   .0000   .6423
ZCOR:    .9338  1.2032  1.0353  1.1862  1.0036
AT% :     .051   1.379    .095    .000  64.150

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
from datetime import datetime


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
relPrj = "/dtsa2Scripts/sim-srm_1155/"
datDir = gitHomeDir + relPrj
jmg.ensureDir(datDir)
print(gitHomeDir)
print(datDir)

srm1555_path = datDir + "sim-srm_1155.msa"


"""
relPrj = "/Documents/work/dtsa2/sim-C-on-Si-7kV"
datDir = homDir + relPrj + "/msa-%g" % (nTraj)
jmg.ensureDir(datDir)
rptDir = homDir + relPrj + "/sim-C-on-Si-7kV Results"
csvFil = homDir + relPrj + "/dtsa2-C-on-Si-%g-kV-kratios-%g-traj.csv" % (e0, nTraj)
"""

start = time.time()
now = datetime.now()
date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
print("start date and time:", date_time) 


DataManager.clearSpectrumList()

"""
St 4591 SRM 1155, Cr-Ni-Mo (AISI 316)
TakeOff = 35.0  KiloVolt = 15.0  Density =  7.970

ELEM:        C      Mn       P       S      Si      Cu      Ni      Cr
XRAY:      ka      ka      ka      ka      ka      ka      ka      ka
ELWT:     .046   1.630    .020    .018    .502    .169  12.180  18.450
KFAC:    .0001   .0163   .0002   .0002   .0033   .0016   .1166   .2128
ZCOR:   4.3792   .9981  1.3232  1.1449  1.5146  1.0843  1.0446   .8669
AT% :     .213   1.649    .036    .031    .993    .148  11.531  19.723

ELEM:        V      Mo      Co      Pb      Fe
XRAY:      ka      la      ka      ma      ka
ELWT:     .047   2.380    .101    .001  64.456
KFAC:    .0005   .0198   .0010   .0000   .6423
ZCOR:    .9338  1.2032  1.0353  1.1862  1.0036
AT% :     .051   1.379    .095    .000  64.150

ELEM:      C              Mn           P                S             Si            Cu            Ni               Cr
ELWT:     0.046   1.630    0.020    0.018  0.502  0 .169  12.180  18.450

ELEM:      V              Mo           Co           Pb           Fe
ELWT:     0.047   2.380    0.101    0.001  64.456


"""


l_comps = [epq.Element.C,  epq.Element.Mn, epq.Element.P,
           epq.Element.S,  epq.Element.Si, epq.Element.Cu,
           epq.Element.Ni, epq.Element.Cr, epq.Element.V,
           epq.Element.Mo, epq.Element.Co, epq.Element.Pb,
           epq.Element.Fe]


l_mfs   = [ 0.046,           1.630,         0.020,
            0.018,           0.502,         0.169,
           12.180,          18.450,         0.047,
            2.380,           0.10,          0.001,
           64.457]

srm1155 = jmg.create_material_from_mf(l_comps, l_mfs, 7.97, "SRM1155")

c  = material("C",  density=2.267)
mn = material("Mn", density=7.43)
p  = material("P",  density=1.88)
s  = material("S",  density=2.0)
si = material("Si", density=2.32)
cu = material("Cu", density=8.96)
mo = material("Mo", density=10.2)
co = material("Co", density=8.86)
pb = material("Pb", density=11.34)
fe = material("Fe", density=7.874)

xrts = []

trs = mc3.suggestTransitions(srm1155, e0)
for tr in trs:
    xrts.append(tr)

srm1155_spc = jm3.simBulkStd(srm1155, det, e0, nTraj, lt, pc, True)
srm1155_spc.setAsStandard(srm1155)
srm1155_spc.display()
srm1155_spc.save(srm1555_path)

"""

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


"""



print "Done!"

date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
print("finish date and time:", date_time) 

end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg




