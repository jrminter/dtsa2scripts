# -*- coding: utf-8 -*-
"""
test-nr-k411-apatite.py

Test N. Ritchie's K411 example
Author: Nicholas W. M. Ritchie 

Scripted by JRM: 2018-10-26
Note: tabs set to 4 spaces

Output with 50,000 trajectories:

K411
(O  All,      0.62222±0.00013),
(Mg All,      0.11069±0.00003),
(Si All,      0.18068±0.00003),
(Ca K-family, 0.9299±0.00034),
(Fe Kα,       0.09731±0.00007),
(Fe Kβ,       0.0977±0.0003),
(Fe L-family, 0.03333±0.00008)

K411 std-20-kV = [O (0.2967±0.0354 mass frac),
                  Mg(0.0784±0.0022 mass frac),
                  Si(0.2400±0.0052 mass frac),
                  Ca(0.1094±0.0011 mass frac),
                  Fe(0.1149±0.0004 mass frac),
                  Σ=0.8394±0.0444]

Note the low total

Done!
This script required 183.110 min
...or 3.052 hr
Elapse: 3:03:06.6

"""

import sys
import os
import time
import shutil
import gov.nist.microanalysis.Utility as epu
import gov.nist.microanalysis.EPQLibrary as epq
import dtsa2.mcSimulate3 as mc3
import dtsa2.jmGen as jmg
import java.util as ju

gitDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/utility"
prjDir = gitDir + relPrj
rptDir = prjDir + '/test-nr-k411-apatite Results/'
spcDir = gitDir + '/dtsa2Scripts/spc/nr-k411-apatite'
jmg.ensureDir(spcDir)


DataManager.clearSpectrumList()
start = time.time()

det   = findDetector("Bruker 5 eV")
e0    =    20.0
nTraj = 50000
dose  =  5000   # nA-sec

def simulate_spectrum(mat, e0, det, nTraj, dose, spcDir):
    """
    simulate_spectrum(mat, e0, det, nTraj, dose, spcDir)
    """
    xrts = []
    trs = mc3.suggestTransitions(mat, e0)
    for tr in trs:
        xrts.append(tr)

    xtraParams={}
    xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))
    # xtraParams.update(mc3.configureOutput(simDir))
    spc = mc3.simulate(mat, det, e0, dose, True, nTraj, True, True, xtraParams)
    sName = "%s std-%g-kV" % (mat.getName(), e0)
    spc.rename(sName)
    spc.setAsStandard(mat)
    spc.display()
    # Save the spectra if NTraj > 500
    if nTraj > 500:
        fi =  spcDir + "/"
        fi += sName
        fi += "-%g-Traj.msa" % (nTraj)
        spc.save(fi)
    return (spc)


mgo = material("MgO", density=3.6)
mgo_spc = simulate_spectrum(mgo, e0, det, nTraj, dose, spcDir)

si = material("Si", density=2.648)
si_spc = simulate_spectrum(si, e0, det, nTraj, dose, spcDir)

apatite = material("Ca(PO4)3", density=3.22)
apatite.setName("Apatite")
apatite_spc = simulate_spectrum(apatite, e0, det, nTraj, dose, spcDir)

fe = material("Fe", density=7.87)
fe_spc = simulate_spectrum(fe, e0, det, nTraj, dose, spcDir)

# this one from database...
k411 = material("K411", density= 2.946)
k411_spc = simulate_spectrum(k411, e0, det, nTraj, dose, spcDir)
print(k411)

res = kratios(k411_spc,{"Si":si_spc,
                        "Mg":mgo_spc,
                        "Ca":apatite_spc,
                        "Fe":fe_spc,
                        "O" :apatite_spc})
print (wrap(res.unknown).kratios())

res_k411 = quantify(k411_spc,{"Si":si_spc,
                              "Mg":mgo_spc,
                              "Ca":apatite_spc,
                              "Fe":fe_spc,
                              "O" :apatite_spc})
# print("dir of qus obj")
# print(dir(res411_1 ))
tmp1 = res_k411.getComposition()
# print("dir of composition obj")
# print(dir(tmp1))
print(tmp1.descriptiveString(False))
k411_spc.setAsMicroanalyticalComposition(tmp1)
k411_spc.display()
if nTraj > 500:
    fi = spcDir + "/K411.xml"
    k411_spc.toXML(fi)


# clean up cruft
shutil.rmtree(rptDir)
print "Done!"

end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg




