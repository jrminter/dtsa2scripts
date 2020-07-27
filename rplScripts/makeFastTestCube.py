# -*- coding: utf-8 -*-
#         1         2         3         4         5         6         7 |
# 23456789012345678901234567890123456789012345678901234567890123456789012
# 
# makeFastTestCube.py
#
# Do (fast) analytical simulations of Cu and Al. Use these to construct
# a data cube with Al in center and Cu on the outside. This is as really
# crude but as fast as it gets. Useful for testing the BTSA-II RPL
# functions. Note how this works in the repo directory.
#
# Date        Ver  Who  Notes
# 2015-05-21 0.90  JRM  Initial example. Verified with Iona v.2015-05-01
#
import dtsa2.hyperTools as ht
import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQTools as et
import math
import os
import shutil


gitDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/rplScripts"
datDir = gitDir + relPrj

os.chdir(datDir)

pyrDir = datDir + "/makeFastTestCube Results"

e0 = 20.0
pc = 1.0
lt = 30.0
det = findDetector("Oxford p4 05eV 2K")
cw   = det.getChannelWidth()
dp   = det.getProperties()
resn = dp.getNumericProperty(epq.SpectrumProperties.Resolution)

rplFil = datDir + "/CuAl.rpl"
rawFil = datDir + "/CuAl.raw"


DataManager.clearSpectrumList()
# make the materials
cu  = epq.Material(epq.Composition([epq.Element.Cu], [1.0]),
                   epq.ToSI.gPerCC(8.096))
cuSpc = simulate(cu, det, keV=e0, dose=lt*pc, withPoisson=True)
cuSpc.rename("ana sim Cu")
al  = epq.Material(epq.Composition([epq.Element.Al], [1.0]),
                   epq.ToSI.gPerCC(2.70))
alSpc = simulate(al, det, keV=e0, dose=lt*pc, withPoisson=True)
alSpc.rename("ana sim Al")

# display(cuSpc)
# display(alSpc)

res=et.RippleFile(64, 64, 2048, et.RippleFile.UNSIGNED , 4, et.RippleFile.BIG_ENDIAN, rplFil, rawFil)
for x in range(-32,32,1):
  print "Working on row %d" % (x)
  for y in range(-32,32,1):
    if terminated:
      break
    r = math.sqrt(x*x+y*y)
    if (r < 16):
      spec = alSpc
    else:
      spec = cuSpc
    cor  = wrap(spec.applyLLT())
    fix  = wrap(cor.positiveDefinite())
    arr = epq.SpectrumUtils.toIntArray(fix)
    # write to Ripple file
    res.write(arr)
    if (x%16==0) and (y%16==0):
      props=cor.getProperties()
      props.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName, "Pixel[%d,%d]" % (x,y))
      props.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
      props.setNumericProperty(epq.SpectrumProperties.FaradayBegin, pc)
      props.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
      display(cor)

# close the rpl files
res.close()      
# append spectrum calibration information to the .rpl file
f=open(rplFil, 'a')
strLine = 'ev-per-chan\t%.4f\n' % (cw)
f.write(strLine)
strLine = 'detector-peak-width-ev\t%.1f\n' % (resn)
f.write(strLine)
f.close()


# clean up cruft
shutil.rmtree(pyrDir)
print "Done!"
