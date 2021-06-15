# sim-Si.py
#   Date      Who  Comment
# ----------  ---  -----------------------------------------------
# 2021-06-14  JRM  sim-Si

import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3","CharacteristicXRayGeneration3,BremsstrahlungXRayGeneration3,FluorescenceXRayGeneration3, XRayTransport3", None)
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
rptDir = wrkDir + '/sim-Si Results/'
spcDir = wrkDir + "/dat/spc"
jmg.ensureDir(spcDir)

os.chdir(wrkDir)

det       = findDetector("Oxford p4 05eV 4K")
e0        =    20    # kV
nTraj     = 10000    # trajectories
lt        =   200    # sec
pc        =     2.0  # nA
imgSzUm   =     5.0  # physical size of images in microns
imgSizePx =   512    # size of images in pixels
vmrlEl    =    40    # number of el for VMRL
dose      = pc * lt  # nA sec

DataManager.clearSpectrumList()

si   = material("Si", density=2.3290)

# Sim Si
xrts = [transition("Si K-L2"), transition("Si K-L3")]
xtraParams={}
xtraParams.update(mc3.configureXRayAccumulators(xrts, charAccum=True, charFluorAccum=True, bremFluorAccum=True))
# note that the image size on the specimen is in meters...
xtraParams.update(mc3.configureEmissionImages(xrts, imgSzUm*1.0e-6, imgSizePx))
xtraParams.update(mc3.configurePhiRhoZ(imgSzUm*1.0e-6))
xtraParams.update(mc3.configureTrajectoryImage(imgSzUm*1.0e-6, imgSizePx))
xtraParams.update(mc3.configureVRML(nElectrons = vmrlEl))
xtraParams.update(mc3.configureOutput(spcDir))

si_spc = mc3.simulate(si, det, e0, dose, nTraj, True, True, xtraParams)
si_spc.rename("Si")
si_spc.setAsStandard(si)
display(si_spc)
fi = spcDir + "/Si-%g-kV-%g-Traj.msa" % (e0, nTraj)
si_spc.save(fi)

print(nTraj)

end  = time.time()

delta = end-start
print(delta)