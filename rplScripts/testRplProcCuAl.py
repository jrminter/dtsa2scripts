# testRplProcCuAl.py
#
# Test processing a rpl file with DTSA-II
# Note : This uses the CuAl .rpl/.raw test cube
# generated with makeFastTestCube
#
# Date        Ver  Who  Notes
# 2015-05-22 0.90  JRM  Initial example. Verified with Iona v.2015-05-01
#                       Longest part seemed to be computing the max px spec
#
import os
import shutil
import dtsa2.hyperTools as ht
import dtsa2 as dt2
import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQTools as ept

gitDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/rplScripts"
datDir = gitDir + relPrj

os.chdir(datDir)

pyrDir = datDir + "/testRplProcCuAl Results"
rplFil = datDir + "/CuAl.rpl"


e0 = 20.0
pc = 1.0
lt = 30.0

det  = findDetector("Oxford p4 05eV 2K")

# Let's start clean...
DataManager.clearSpectrumList()

rs = ht.openRipple(rplFil, e0, pc, lt, det)
rs.setPosition(32, 32) # go to a position
dt2.display(rs)
ary = wrap(rs).toDouble()
print("Print the min and max")
print(min(ary))
print(max(ary))

dir(rs)
print(dir(rs))

print("Computing maximum pixel spectrum")

mps = wrap(ht.maxPixel(rs))
mps.rename("CuAl max px")
dt2.display(mps)


# clean up cruft
shutil.rmtree(pyrDir)
print "Done!"
