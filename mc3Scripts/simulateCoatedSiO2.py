# simulateCoatedSiO2.py
#   Date      Who  Comment
# ----------  ---  -----------------------------------------------
# 2016-09-12  JRM  Test simulation of a coated SiO2 for TMK


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

rptDir = wrkDir + '/simulateCoatedSiO2 Results/'

spcDir = wrkDir + "/dat/spc"
jmg.ensureDir(spcDir)

det     = findDetector("Oxford p4 05eV 2K")
e0      = 5       # kV
nTraj   = 10000   # trajectories
lt      = 60      # sec
pc      = 0.25    # nA
tNmTa   = 2
tNmC    = 10
tNmAuPd = 10
tNmSiO2 = 100


dose = pc * lt  # na-sec"

DataManager.clearSpectrumList()

sio2 = material("SiO2", density=2.65)
c = material("C", density=2.1)
aupd = epq.Material(epq.Composition([epq.Element.Au, epq.Element.Pd],
                                    [0.60,0.40]),
                                    epq.ToSI.gPerCC(0.6*19.30+0.4*11.9))
ta = material("Ta", density=16.4)

layers = [ [sio2, tNmSiO2*1.0e-9],
           [c, 50.0e-6]
         ]

# multiFilm(layers, det, e0=20.0, withPoisson=True,
#           nTraj=defaultNumTraj, dose=defaultDose,
#           sf=defaultCharFluor, bf=defaultBremFluor,
#            xtraParams=defaultXtraParams)



basSim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                       dose=dose, sf=True, bf=True, xtraParams={})
sName = "%g nm SiO2 on C" % tNmSiO2
basSim.rename(sName)
basSim.display()

fi = spcDir + "/%g-nm-SiO2-on-Si-%g-kV-%g-Traj.msa" % (tNmSiO2, e0, nTraj)
basSim.save(fi)

layers = [ [aupd, tNmAuPd*1.0e-9],
           [sio2, tNmSiO2*1.0e-9],
           [c, 50.0e-6]
         ]


aupdCtd = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                       dose=dose, sf=True, bf=True, xtraParams={})
sName = "%g nm AuPd on %g nm SiO2 on C" % (tNmAuPd, tNmSiO2)
aupdCtd.rename(sName)
aupdCtd.display()
fi = spcDir + "/%g-nm-AuPd-%g-nm-SiO2-on-C-%g-kV-%g-Traj.msa" % (tNmAuPd, tNmSiO2, e0, nTraj)
aupdCtd.save(fi)


layers = [ [c, tNmC*1.0e-9],
           [sio2, tNmSiO2*1.0e-9],
           [c, 50.0e-6]
         ]
cCtd =  mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                       dose=dose, sf=True, bf=True, xtraParams={})
sName = "%g nm C on %g nm SiO2 on C" % (tNmAuPd, tNmSiO2)
cCtd.rename(sName)
cCtd.display()
fi = spcDir + "/%g-nm-C-%g-nm-SiO2-on-C-%g-kV-%g-Traj.msa" % (tNmC, tNmSiO2, e0, nTraj)
cCtd.save(fi)


layers = [ [ta, tNmTa*1.0e-9],
           [sio2, tNmSiO2*1.0e-9],
           [c, 50.0e-6]
         ]

taCtd =  mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                       dose=dose, sf=True, bf=True, xtraParams={})
sName = "%g nm Ta on %g nm SiO2 on C" % (tNmTa, tNmSiO2)
taCtd.rename(sName)
taCtd.display()
fi = spcDir + "/%g-nm-Ta-%g-nm-SiO2-on-C-%g-kV-%g-Traj.msa" % (tNmTa, tNmSiO2, e0, nTraj)
cCtd.save(fi)



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


