"""
testWindowTransmission.py

"""

import os, shutil

homDir  = os.environ['HOME']
homDir  = homDir.replace('\\','/')
relPrj  = "/Documents/git/dtsa2Scripts/utility"
prjDir  = homDir + relPrj
rptDir  = prjDir + '/testWindowTransmission Results/'


det = findDetector("Oxford p4 05eV 2K")
na = det.getName()
na += " Window Transmission"
print(na)

wt = windowTransmission(det)
wt.rename(na)
wt.display()




# clean up cruft
shutil.rmtree(rptDir)
print "Done!"

