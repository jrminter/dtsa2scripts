# simulateCoatedEagleXG.py
#   Date      Who  Comment
# ----------  ---  -----------------------------------------------
# 2017-06-20  JRM  Simulate Corning EagleXG from patent midpoints
#                  Elapse: 4:11:55.1 on crunch for 20000 traj
# 2017-06-21  JRM  Let's try this at 7 kV. F(chi) of B is better
#                  and it should be faster too! Sorted mass fraction
#                  list too. Elapse: 1:16:15.3 on ROC...

import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3,FluorescenceXRayGeneration3, XRayTransport3", None)


import os

import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQLibrary.Detector as epd
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.microanalysis.NISTMonte.Gen3 as nm3
import gov.nist.microanalysis.EPQTools as et
import gov.nist.microanalysis.dtsa2 as dt2

import dtsa2.jmGen as jmg
import dtsa2.mcSimulate3 as mc3
import dtsa2.hyperTools as ht
import dtsa2.jmMC3 as jm3

import shutil
import time

import java.util as jutil
import java.io as jio
import java.nio.charset as cs
import string


start = time.time()

gitHom = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts"
prjDir = gitHom + relPrj
wrkDir = prjDir + "/mc3Scripts"
csvDir = wrkDir + "/dat/csv"
jmg.ensureDir(wrkDir)

os.chdir(wrkDir)

rptDir = wrkDir + '/simulateCoatedEagleXG Results/'

spcDir = wrkDir + "/dat/spc"
simDir = wrkDir + "/dat/sim"
jmg.ensureDir(spcDir)
jmg.ensureDir(simDir)
jmg.ensureDir(spcDir)

print(simDir)

det      = findDetector("Oxford p4 05eV 2K")
e0       =     7     # kV
nTraj    = 20000     # trajectories
lt       =   100     # sec
pc       =     1.0   # nA
tNmC     =    20     # nm C


dose = pc * lt  # na-sec"

DataManager.clearSpectrumList()

eagleXG = epq.Material(epq.Composition([epq.Element.O,
                                        epq.Element.Si,
                                        epq.Element.Al,
                                        epq.Element.Ca,
                                        epq.Element.B,
                                        epq.Element.Mg,
                                        epq.Element.Sr,
                                        epq.Element.Sn,
                                        epq.Element.Ba,
                                        epq.Element.As,
                                        epq.Element.Sb,
                                        epq.Element.Fe,
                                        epq.Element.Zr,
                                        epq.Element.Ti],
                                        [0.518777,
                                         0.301391,
                                         0.090081,
                                         0.038763,
                                         0.032655,
                                         0.007728,
                                         0.006965,
                                         0.001198,
                                         0.001092,
                                         0.000238,
                                         0.000387,
                                         0.000355,
                                         0.000218,
                                         0.000152]),
                      epq.ToSI.gPerCC(2.38))
eagleXG.setName("EagleXG")

a = jm3.simCarbonCoatedStd(eagleXG, det, e0, nTraj, lt, pc, tNmC)
a.display()


# clean up cruft
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




