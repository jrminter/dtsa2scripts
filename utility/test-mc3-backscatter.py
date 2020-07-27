"""
test-mc3-backscatter.py
"""

import os, shutil
import dtsa2.mcSimulate3 as mc3

start = time.time()

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
outDir = homDir + "/Documents/git/dtsa2Scripts/utility/output/"
rptDir = wrkDir + '/test-mc3-backscatter Results/'

e0 = 20.0
zno = material("ZnO", density=5.61)

fName = outDir + "bse-test.csv"

out = mc3.backscatter(zno, e0, 10000)
print(out)


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
