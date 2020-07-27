"""
test-mc3-energy-deposition.py
"""

import os, shutil
import dtsa2.mcSimulate3 as mc3

start = time.time()

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
outDir = homDir + "/Documents/git/dtsa2Scripts/utility/output/"
rptDir = wrkDir + '/test-mc3-energy-deposition Results/'

c   = material("C", density=2.267)
zno = material("ZnO", density=5.61)

layers = [ [c, 20*1.0e-9], [zno, 50*1.0e-6] ]
mc3.energyDeposition(layers, 15, 1000, outDir)


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
