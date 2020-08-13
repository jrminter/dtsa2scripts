# simulateCoatedSiO2.py
#   Date      Who  Comment
# ----------  ---  -----------------------------------------------
# 2020-08-11  JRM  Simulate SiO2 on Si using functions from jmMc3



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
import dtsa2.jmMC3 as jm3
import dtsa2.mcSimulate3 as mc3
import dtsa2.hyperTools as ht

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
bClean = True
jmg.ensureDir(wrkDir)

os.chdir(wrkDir)

rptDir = wrkDir + '/simulateCoatedSiO2 Results/'

spcDir = wrkDir + "/dat/spc"
jmg.ensureDir(spcDir)

det     = findDetector("Oxford p4 05eV 2K")
e0      = 5       # kV
nTraj   = 10000   # trajectories
lt      = 60      # sec
pc      = 0.25    # nA
tNmSiO2 = 20

si = material("Si", 2.329)
sio2 = material("SiO2", 2.65)

DataManager.clearSpectrumList()
#                                                      em img          emi size                                        
siSpc = jm3.uncoatedSimBulkStd(si,det,e0,nTraj,spcDir, 5.0e-6, lt, pc, 512)
siSpc.display()
fi = spcDir + "/Si-std-%g-kV-%g-Traj.msa" % (e0, nTraj)
siSpc.save(fi)

sio2Spc = jm3.uncoatedSimBulkStd(sio2,det,e0,nTraj,spcDir, 5.0e-6, lt, pc, 512)
sio2Spc.display()
fi = spcDir + "/SiO-std-%g-kV-%g-Traj.msa" % (e0, nTraj)
sio2Spc.save(fi)



print "Done!"

end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg


