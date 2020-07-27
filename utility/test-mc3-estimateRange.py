"""
test-mc3-estimateRange.py
"""


import os, shutil
import dtsa2.mcSimulate3 as mc3

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
outDir = homDir + "/Documents/git/dtsa2Scripts/utility/output/"
rptDir = wrkDir + '/test-mc3-estimateRange Results/'


e0 = 20.0
mat = material("Fe2O3",5.0)
res = mc3.estimateRange(mat, e0)
print(res)

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
