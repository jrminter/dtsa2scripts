"""
test-mc3-simulate-coated-block.py
"""


import os, shutil
import dtsa2.mcSimulate3 as mc3

start = time.time()

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
outDir = homDir + "/Documents/git/dtsa2Scripts/utility/output/"
rptDir = wrkDir + '/test-mc3-simulate-coated-block Results/'

det = findDetector("Probe")
e0  = 20.0
ht  = 0.5e-06
wd  = 2.0e-6
ctT = 0.5e-6
cov = 0.5 # coverage
det = findDetector("Probe")
# al2o3 = material("Al2O3", 3.95)
al = material("Al", 2.7)
si = material("Si", 2.328)
zno = material("ZnO", 5.61)


DataManager.clearSpectrumList()

spcCtdBlock = mc3.coatedBlock(zno, ht, wd, al, ctT, si, det, e0,
                              True, 1000, 120.0, True, True, {})


# sName = "%g um Al2O3 spheres on Fe2O3 at %g kV" % (rad*1.0e6, e0)
# spc.rename(sName)
spcCtdBlock.display()

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
