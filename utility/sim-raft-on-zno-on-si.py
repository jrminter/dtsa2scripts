# -*- coding: utf-8 -*-
"""
sim-raft-on-zno-on-si.py

Done!
This script required 14.898 min
Elapse: 0:14:53.9

"""

import sys
import os
import time
import shutil
import gov.nist.microanalysis.Utility as nu
import dtsa2.mcSimulate3 as mc3
import dtsa2.jmGen as jmg
import java.util as ju

# adjustable parameters

det      = findDetector("Oxford p4 05eV 2K")
nTraj    = 50000     # trajectories
e0       =     7.0
tNmRAFT  =    10
tNmZnO   =    50
vmrlEl   =    50
przDepUm =     1.0
imgSzUm  =     1.0
imgSize  =   512    # pix


dose    =  5000   # nA-sec
nDigits = 6


# directories
gitDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/utility"
prjDir = gitDir + relPrj
resDir = prjDir + '/sim-raft-on-zno-on-si'
jmg.ensureDir(resDir)
rptDir = prjDir + '/sim-raft-on-zno-on-si Results/'
spcDir = gitDir + '/dtsa2Scripts/spc/sim-raft-on-zno-on-si'
jmg.ensureDir(spcDir)


DataManager.clearSpectrumList()
start = time.time()

# define materials
zno  = material("ZnO", density=5.61)
si   = material("Si", density=2.33)
raft = material("C17H35O5S2P1", density=2.0)
raft.setName("RAFT polymer")

# specify the x-ray transitions
xrts = []
trs = mc3.suggestTransitions(zno, e0)
for tr in trs:
    xrts.append(tr)
trs = mc3.suggestTransitions(si, e0)
for tr in trs:
    xrts.append(tr)
trs = mc3.suggestTransitions(raft, e0)
for tr in trs:
    xrts.append(tr)

xtraParams={}
xtraParams.update(mc3.configurePhiRhoZ(przDepUm*1.0e-6))
xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))
# note that the image size on the specimen is in meters...
xtraParams.update(mc3.configureEmissionImages(xrts, imgSzUm*1.0e-6, imgSize))
xtraParams.update(mc3.configureTrajectoryImage(imgSzUm*1.0e-6, imgSize))
xtraParams.update(mc3.configureVRML(nElectrons = vmrlEl))
xtraParams.update(mc3.configureOutput(spcDir))
print(xtraParams)

layers = [ [raft, tNmRAFT*1.0e-9], # top layer of RAFT
           [zno,   tNmZnO*1.0e-9], # ZnO layer
           [raft, tNmRAFT*1.0e-9], # bottom layer of RAFT
           [si,   50.0e-6]         # 50 micron (inf. thick Si substrate)
         ]

spc = mc3.multiFilm(layers, det, e0, True, nTraj,
                      dose, True, True, xtraParams)
strName = "%g-nm-RAFT-encapsulated-%g-nm-ZnO-on-Si-%g-kV"
spcName =  strName % (tNmRAFT, tNmZnO, e0)

spc.rename(spcName)
spc.display()

# only save spectrum with > 500 traj
if nTraj > 500:
    fi =  spcDir + "/"
    fi += spcName
    fi += "-%g-Traj.msa" % (nTraj)
    spc.save(fi)






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



