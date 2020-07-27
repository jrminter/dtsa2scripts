# testRplProcPaint.py
#
# Test processing a rpl file with DTSA-II
# Note : rpl is converted from Oxford Paint set:
# 1. To Inca project, exported to .rpl/.raw Ripple file pair.
# 2. Converted with hyperspy to vector.
# 3. .rpl cleaned up manually... (need to fix...)
#
# Create the Oxford-Paint detector (just add it to the default probe
# instrument) by importing the oxford.msa spectrum from DTSA's add
# detector dialog. You need to manually set the size to 80 mm2 and
# the window to Moxtek manufacturers.
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

pyrDir = datDir + "/testRplProcPaint Results"
rplFil = datDir + "/paint-vec.rpl"


e0 = 20.0
pc = 1.0
lt = 30.0

det  = findDetector("Oxford-Paint")

if(det.toString()[0:4] != 'Oxfo'):
	print("You need to add the Oxford-Paint detector. It is using the default...")

# Let's start clean...
DataManager.clearSpectrumList()

rs = ht.openRipple(rplFil, e0, pc, lt, det)
rs.setPosition(111, 35) # go to a position
dt2.display(rs)
ary = wrap(rs).toDouble()
print("Print the min and max")
print(min(ary))
print(max(ary))

dir(rs)
print(dir(rs))

print("Computing maximum pixel spectrum")

mps = wrap(ht.maxPixel(rs))
mps.rename("Paint maximum px")
dt2.display(mps)


# clean up cruft
shutil.rmtree(pyrDir)
print "Done!"
