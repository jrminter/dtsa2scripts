"""
test-mc3-simulate-coated-substrate.py
"""


import os, shutil
import dtsa2.mcSimulate3 as mc3

start = time.time()

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
outDir = homDir + "/Documents/git/dtsa2Scripts/utility/output/"
rptDir = wrkDir + '/test-mc3-simulate-coated-substrate Results/'

e0  = 20.0
det = findDetector("Probe")
mat = material("Fe2O3",5.0)
ctg = material("C", 1.8)

DataManager.clearSpectrumList()

spc = mc3.coatedSubstrate(ctg, 0.020, mat, det, 20.0,
                          True, 1000, 120.0, True,
                          True, {})
sName = "20 nm C on Fe2O3 at 20kV"
spc.rename(sName)
spc.display()

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
