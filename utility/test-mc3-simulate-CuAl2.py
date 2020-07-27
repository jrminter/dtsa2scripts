"""
test-mc3-simulate-CuAl2.py
"""


import os, shutil
import dtsa2.mcSimulate3 as mc3

start = time.time()

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
outDir = homDir + "/Documents/git/dtsa2Scripts/utility/output/"
rptDir = wrkDir + '/test-mc3-simulate-CuAl2 Results/'

DataManager.clearSpectrumList()

e0  = 15.0
det = findDetector("Probe")
mat = material("CuAl2", 4.36)

spc = mc3.simulate(mat, det, e0, 120.0, True, 10000, True, True, {})
sName = "Bulk CuAl2"
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
