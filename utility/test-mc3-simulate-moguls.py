"""
test-mc3-simulate-moguls.py
"""


import os, shutil
import dtsa2.mcSimulate3 as mc3

start = time.time()

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
outDir = homDir + "/Documents/git/dtsa2Scripts/utility/output/"
rptDir = wrkDir + '/test-mc3-simulate-moguls Results/'

det = findDetector("Probe")
e0  = 20.0
rad = 5.0e-06 # micron sphere
ht  = 0.5 * rad
cov = 0.5 # coverage
det = findDetector("Probe")
mat = material("Fe2O3", 3.95)


DataManager.clearSpectrumList()

spcMoguls = mc3.moguls(mat, ht, True, cov, det, e0, True, 1000, 120.0, True,
                       True, {})


# sName = "%g um Al2O3 spheres on Fe2O3 at %g kV" % (rad*1.0e6, e0)
# spc.rename(sName)
spcMoguls.display()

spcBulk = mc3.simulate(mat, det, e0, 120.0, True, 1000, True, True, {})
spcBulk.display()

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
