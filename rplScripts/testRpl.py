import os
import sys
import shutil
import fnmatch
import math
import dtsa2 as dt2
import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQTools as ept
from java.lang import Double

import dtsa2.jmGen as jmg


gitRoot = os.environ['GIT_HOME']
datDir  = gitRoot + "/dtsa2Scripts/rplScripts"
pyrDir="./testRpl Results"
sumSpec = datDir + "/paint-sum-dtsa.msa"
rplFil = datDir + "/paint-vec.rpl"
rawFil = datDir + "/paint-vec.raw"
os.chdir(datDir)

print(os.getcwd())
DataManager.clearSpectrumList()
# calibrated with extracted BaSO4 spectrum
myDet = findDetector("Oxford-Paint")

dc = ept.RippleFile(rplFil, True)

# 111, 35 is the location of the center of the BaSO4 particle

mySpc = jmg.getSpectrumFromDataCube(dc, myDet, sumSpec, 111, 35, 14387., pc=1.0, bDebug=False)


dc.close() 

# clean up cruft
shutil.rmtree(pyrDir)
print "Done!"
