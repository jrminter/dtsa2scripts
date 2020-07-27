"""
test-mc3-simulate-sphere-bed.py
"""


import os, shutil
import dtsa2.mcSimulate3 as mc3

start = time.time()

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
outDir = homDir + "/Documents/git/dtsa2Scripts/utility/output/"
rptDir = wrkDir + '/test-mc3-simulate-sphere-bed Results/'

e0  = 20.0
rad = 0.5e-06 # micron sphere
det = findDetector("Probe")
sub = material("Fe2O3",5.0)
sph = material("Al2O3", 3.95)


DataManager.clearSpectrumList()

spc = mc3.sphereBed(sph, rad, sub, False, det,
                    e0, True, 1000, 120.0, True, True, {})
sName = "%g um Al2O3 spheres on Fe2O3 at %g kV" % (rad*1.0e6, e0)
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
