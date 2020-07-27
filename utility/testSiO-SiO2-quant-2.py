# This Python file uses the following encoding: utf-8
# testSiO-SiO2-quant-2.py
# This version reads from disk
# 2018-10-16

import sys
import os
import time
import shutil

# import dtsa2 as dt2
import dtsa2.mcSimulate3 as mc3
import dtsa2.jmGen as jmg

gitDir    = os.environ['GIT_HOME']
relPrj    = "/dtsa2Scripts/utility"
prjDir    = gitDir + relPrj
rptDir    = prjDir + '/testSiO-SiO2-quant-2 Results/'
spcDir    = gitDir + relPrj + "/sim-quant-sio-w-sio2"

det       = findDetector("Oxford p4 05eV 2K")
e0        =     4.0   # kV
nDigits   =     5

DataManager.clearSpectrumList()

start = time.time()

sio2 = material("SiO2", density=2.196)
# sio  = material("SiO", density=2.13)

# Read standard
fi = spcDir + "/SiO2 std.msa"
spc_sio2_std = readSpectrum(fi)
spc_sio2_std.display()

lDoseUnk = [50, 100, 200, 500, 1000, 2500, 5000]
xrts = []

trs = mc3.suggestTransitions(sio2, e0)
for tr in trs:
    xrts.append(tr)

stds =  { element("O"): spc_sio2_std, element("Si"): spc_sio2_std }

qus = multiQuant(det, e0, stds, {})

for doseUnk in lDoseUnk:
    sName = "SiO Unk %g nA-sec.msa" % (doseUnk)
    fi = spcDir + "/" + sName
    spc_sio_unk = readSpectrum(fi)
    spc_sio_unk.display()
    res = qus.compute(spc_sio_unk)
    print(sName)
    print("Weight Fraction")
    print("Si")
    siWF = res.getComposition().weightFractionU(element("Si"), True)
    print(jmg.pretty_print_unc_val(siWF, nDigits))
    print("O")
    oWF = res.getComposition().weightFractionU(element("O"), True)
    print(jmg.pretty_print_unc_val(oWF, nDigits))
    print("")
    print("Atomic Percent")
    print("Si")
    siAP = res.getComposition().atomicPercentU(element("Si"), True)
    print(jmg.pretty_print_unc_val(siAP, nDigits))
    print("O")
    oAP = res.getComposition().atomicPercentU(element("O"), True)
    print(jmg.pretty_print_unc_val(oAP, nDigits))

    print("\n")


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