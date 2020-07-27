"""
test-mc3-zaf.py

"""

import os, shutil
import dtsa2.mcSimulate3 as mc3

start = time.time()

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
rptDir = wrkDir + '/test-mc3-zaf Results/'

mc3.zaf(material("Al2O3",1),d2,10.0,stds={ "Al":"Al", "O":"MgO" })


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
