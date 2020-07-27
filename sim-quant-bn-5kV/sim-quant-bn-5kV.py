# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

sim-BN-quant.py

A script to model and quantify BN at 5 kV

  Date      Who   Ver    Comment
----------  ---  -----   ------------------------------------------
2020-05-24  JRM  0.0.1   Initial test. There is an issue with 
                          __init.py__
                          10,000 trajectories
                          BN-5kV = [B(0.4575±0.0963 mass frac)
                                    N(0.5825±0.0883 mass frac),
                                    Σ=1.0399±0.1846]
Done!
This script required 5.020 min
"""

__revision__ = "$Id: sim-quant-bn-5kV.py John R. Minter $"
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

det = findDetector("Oxford p4 05eV 2K")
e0       =     5    # kV
nTraj    =  20000   # trajectories
lt       =    500   # sec
pc       =    5.0   # nA
dose     = pc * lt  # na-sec"
bSaveSpc = True
coating = False

datDir = ""


gitHomeDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/sim-quant-bn-5kV"
datDir = gitHomeDir + relPrj
print(gitHomeDir)
print(datDir)

start = time.time()

DataManager.clearSpectrumList()

# define materials
b  = material("B", density = 2.37)
bn = material("BN", density = 2.1)
aln = material("AlN", density = 3.26)

xrts = []

trs = mc3.suggestTransitions(bn, e0)
for tr in trs:
  xrts.append(tr)

# def simulateBulkStandard(mat, name, det, e0, lt, pc, withPoisson=True, nTraj=100, sf=True, bf=True, xtraParams={}):

spc_b = jm3.simBulkStd(b, det, e0, nTraj, 100, 1.0, False)
spc_b.display()
spc_b.rename("B-5kV")
spc_b.setAsStandard(b)
fi = datDir + "/spc_b.msa"
print(fi)
if(bSaveSpc):
  spc_b.save(fi)

spc_bn = jm3.simBulkStd(bn, det, e0, nTraj, 100, 1.0, False)
spc_bn.display()
spc_bn.rename("BN-5kV")
fi = datDir + "/spc_bn.msa"
print(fi)
if(bSaveSpc):
  spc_bn.save(fi)

spc_aln = jm3.simBulkStd(aln, det, e0, nTraj, 100, 1.0, False)
spc_aln.display()
spc_aln.rename("AlN-5kV")
spc_aln.setAsStandard(aln)
fi = datDir + "/spc_aln.msa"
print(fi)
if(bSaveSpc):
  spc_aln.save(fi)

qua = quantify(spc_bn, {"B":spc_b, "N":spc_aln},
               refs={}, preferred=(), elmByDiff=None,
               oByStoic=False, oxidizer=None, extraKRatios=None, fiat={})
tmp = qua.getComposition()
print(tmp.descriptiveString(False))


print("%d trajectory analysis Done!" % nTraj)

end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg


