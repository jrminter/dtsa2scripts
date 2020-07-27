# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

sim-BN-quant.py

A reproducible example to model and quantify BN at 5 kV

  Date      Who   Ver    Comment
----------  ---  -----   ------------------------------------------
2020-05-24  JRM  0.0.1   Initial test. There is an issue with 
                         Note: This fails...
"""


__revision__ = "$Id: sim-quant-bn-5kV.py John R. Minter $"
__version__ = "0.0.1"


import dtsa2 as dt2
import gov.nist.microanalysis.dtsa2 as gdtsa2
import dtsa2.mcSimulate3 as mc3

# change the datDir to match your system
# MacOS
# datDir   = "/Users/jrminter/dat/results"
# Windows10
datDir =  "D:/dtsa2_sims/sim-bn-quant"

det = findDetector("Bruker 5 eV")
e0       =     2.5  # kV
nTraj    =  10000   # trajectories
lt       =   500    # sec
pc       =     5.0  # nA
dose     = pc * lt  # na-sec"
bSaveSpc = True
# coating  = False


def simBulkStd(mat, det, e0, nTraj, lt=100, pc=1.0, ctd=True):
    """simBulkStd(mat, det, e0, nTraj, lt=100, pc=1.0, ctd=True)

    Use mc3 simulation to simulate an uncoated standard specimen

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
    ctd - Boolean (True) - is C coated


    Returns
    -------
    sim - DTSA scriptable spectrum 
        The simulated standard spectrum
    
    Example
    -------
    det = findDetector("Brucker 5eV")
    cu = material("Cu", density=8.92)
    a = simBulkStd(cu, det, 20.0, 100, 100, 1.0)
    a.display()
    """
    
    dose = pc * lt  # na-sec"
    xrts = []

    trs = mc3.suggestTransitions(mat, e0)
    for tr in trs:
        xrts.append(tr)

    xtraParams={}
    xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))

    sim = mc3.simulate(mat, det, e0, dose, withPoisson=True, nTraj=nTraj,
                       sf=True, bf=True, xtraParams=xtraParams)

    sName = "%s-%g-kV" % (mat, e0)
    sim.rename(sName)
    sim.setAsStandard(mat)
    return sim

start = time.time()

DataManager.clearSpectrumList()

# define materials
b  = material("B", density = 2.37)
bn = material("BN", density = 2.1)
aln = material("AlN", density = 3.26)

xrts = []

trs = mc3.suggestTransitions(b, e0)
for tr in trs:
  xrts.append(tr)
trs = mc3.suggestTransitions(aln, e0)
for tr in trs:
  xrts.append(tr)

spc_b = simBulkStd(b, det, e0, nTraj, 100, 1.0, False)
spc_b.display()
spc_b.rename("B-5kV")
spc_b.setAsStandard(b)
fi = datDir + "/spc_b.msa"
print(fi)
if(bSaveSpc):
  spc_b.save(fi)

spc_aln = simBulkStd(aln, det, e0, nTraj, 100, 1.0, False)
spc_aln.rename("AlN-5kV")
spc_aln.setAsStandard(aln)
spc_aln.display()
fi = datDir + "/spc_aln.msa"
print(fi)
if(bSaveSpc):
  spc_aln.save(fi)

# Sim BN "unknown"
spc_bn = simBulkStd(bn, det, e0, nTraj, 100, 1.0, False)
spc_bn.rename("BN-5kV")
spc_bn.setAsStandard(bn)
spc_bn.display()
fi = datDir + "/spc_bn.msa"
print(fi)
if(bSaveSpc):
  spc_bn.save(fi)

# Run a quantification - that fails...

qua = quantify(spc_bn, {"B":spc_b, "N":spc_aln},
               refs={}, preferred=(), elmByDiff=None,
               oByStoic=False, oxidizer=None, extraKRatios=None, fiat={})
tmp = qua.getComposition()
print(tmp.descriptiveString(False))

end = time.time()

delta = end - start

print(delta)



