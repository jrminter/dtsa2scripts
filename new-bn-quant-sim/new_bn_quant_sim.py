# -*- coding: utf-8 -*-

"""
/Users/jrminter/dat/dtsa2-eds-sim/sim-bn-quant/new_bn_quant_sim.py

Runs with Kelvin 2018-09-26
Running /Users/jrminter/dat/dtsa2-eds-sim/sim-bn-quant/new_bn_quant_sim.py
BN-5-kV = [B(0.4968±0.0957 mass frac),N(0.5944±0.0850 mass frac),Σ=1.0912±0.1807]

Runs with Lorentz 2019-07-29
Running /Users/jrminter/dat/dtsa2-eds-sim/sim-bn-quant/new_bn_quant_sim.py
BN-5-kV = [B(0.4887±0.0932 mass frac),N(0.5849±0.0831 mass frac),Σ=1.0736±0.1762]

Crashes with latest Lorentz with
File "/Applications/NIST DTSA-II Lorentz 2020-05-26/Lib/dtsa2/__init__.py", line 1320, in quantify
    return qus.compute(unknown)
  at Jama.SingularValueDecomposition.<init>(SingularValueDecomposition.java:181)



"""

__revision__ = "$Id: sim-quant-bn.py John R. Minter $"
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
import dtsa2.mcSimulate3 as mc3

det = findDetector("Bruker 5 eV")
e0       =     5    # kV
nTraj    =  10000   # trajectories
lt       =    500   # sec
pc       =    5.0   # nA
dose     = pc * lt  # na-sec"
bSaveSpc = True

def simBulkStandard(mat, name, det, e0, lt, pc, nTraj,
                    withPoisson=True, sf=True, bf=True, xtraParams={}):
    """simBulkStandard(mat, name, det, e0, lt, pc, nTraj,
    withPoisson=True, sf=True, bf=True, xtraParams={})"""

    std = mc3.simulate(mat, det, e0, lt*pc, withPoisson=True,
                       nTraj=nTraj, sf=True, bf=True, xtraParams={})
    props=std.getProperties()
    props.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
    props.setNumericProperty(epq.SpectrumProperties.FaradayBegin, pc)
    props.setNumericProperty(epq.SpectrumProperties.FaradayEnd, pc)
    props.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
    std.setAsStandard(mat)
    return(std)

start = time.time()

DataManager.clearSpectrumList()

# define materials
b  = material("B", density = 2.37)
bn = material("BN", density = 2.1)
aln = material("AlN", density = 3.26)


b_spec = simBulkStandard(b, "B-std", det, e0, lt, pc, nTraj,
                         withPoisson=True, sf=True, bf=True,
                         xtraParams={})
b_spec.display()
b_spec.rename("B-%g-kV" % e0)

aln_spec = simBulkStandard(aln, "AlN-std", det, e0, lt, pc, nTraj,
                           withPoisson=True, sf=True, bf=True,
                           xtraParams={})
aln_spec.display()
aln_spec.rename("AlN-%g-kV" % e0)

# This does not help either...

bn_spec = simBulkStandard(bn, "BN-std", det, e0, lt, pc, nTraj,
                          withPoisson=True, sf=True, bf=True,
                          xtraParams={})
bn_spec.display()
bn_spec.rename("BN-%g-kV" % e0)

xrts = []

trs = mc3.suggestTransitions(b, e0)
for tr in trs:
  xrts.append(tr)

trs = mc3.suggestTransitions(aln, e0)
for tr in trs:
  xrts.append(tr)

qua = quantify(bn_spec, {"B":b_spec, "N":aln_spec})
tmp = qua.getComposition()
print(tmp.descriptiveString(False))

# result with std bundles
# BN-5-kV = [B(0.5014±0.0972 mass frac),N(0.6001±0.0861 mass frac),Σ=1.1015±0.1833]

print "Done!"






end = time.time()

elapsed = end - start
print(elapsed)
