"""
test-jmc3-simulate-bulk.py

Try to extract emitted and detected spectra
"""


import os, shutil
import dtsa2.jmcSimulate3 as jmc3

start = time.time()

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
outDir = homDir + "/Documents/git/dtsa2Scripts/utility/output/"
rptDir = wrkDir + '/test-jmc3-simulate-bulk Results/'

DataManager.clearSpectrumList()

e0  = 20.0
det = findDetector("Probe")
mat = material("Fe2O3",5.0)

lSpec = jmc3.simulate(mat, det, e0, 120.0, True, 100, True, True, {})
sName = "Bulk Fe2O3 detected"
lSpec[0].rename(sName)
lSpec[0].display()
sName = "Bulk Fe2O3 emitted"
lSpec[1].rename(sName)
lSpec[1].display()

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
