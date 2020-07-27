# testSiO-SiO2-quant.py

import sys
import os
import time
import shutil

# import dtsa2 as dt2
import dtsa2.mcSimulate3 as mc3

gitDir    = os.environ['GIT_HOME']
relPrj    = "/dtsa2Scripts/utility"
prjDir    = gitDir + relPrj
rptDir    = prjDir + '/testSiO-SiO2-quant Results/'

det       = findDetector("Oxford p4 05eV 2K")
e0        =     4.0   # kV
nTrajStd  = 50000     # trajectories
nTrajUnk  = 50000     # trajectories
doseStd   =  5000     # nA-sec

DataManager.clearSpectrumList()

start = time.time()

sio2 = material("SiO2", density=2.196)
sio  = material("SiO", density=2.13)

lDoseUnk = [50, 100, 200, 500, 1000, 2500, 5000]


xrts = []

trs = mc3.suggestTransitions(sio2, e0)
for tr in trs:
    xrts.append(tr)



xtraParams={}
xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))
# xtraParams.update(mc3.configureOutput(simDir))

print(xtraParams)

spc_sio2_std = mc3.simulate(sio2, det, e0, doseStd, True, nTrajStd, True, True, xtraParams)
sName = "SiO2 std"
spc_sio2_std.rename(sName)
spc_sio2_std.setAsStandard(sio2)
spc_sio2_std.display()

stds =  { element("O"): spc_sio2_std, element("Si"): spc_sio2_std }

# oStd  = {"El":element("O"),  "Spc":spc_sio2_std}
# 
# stds = [oStd, siStd]

for doseUnk in lDoseUnk:
    spc_sio = mc3.simulate(sio, det, e0, doseUnk, True, nTrajStd, True, True, xtraParams)
    sName = "SiO Unk %g nA-sec" % (doseUnk)
    spc_sio.rename(sName)
    spc_sio.display()
    res = quantify(spc_sio, stds, {}, (), None, False, None, None, {})
    print(res.getComposition().toHTMLTable())


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