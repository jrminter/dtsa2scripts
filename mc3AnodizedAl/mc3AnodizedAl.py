# -*- coding: utf-8 -*-
# DTSA-II Script - J. R. Minter - 2016-04-18

import sys
import os
import glob
import shutil
import time
import math
import csv

sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3, FluorescenceXRayGeneration3, XRayTransport3", None)
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.microanalysis.NISTMonte.Gen3 as nm3
import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQLibrary.Detector as epd
import gov.nist.microanalysis.Utility as epu
import gov.nist.microanalysis.EPQTools as ept


import dtsa2 as dt2
import dtsa2.mcSimulate3 as mc3
import dtsa2.jmGen as jmg

"""A series of wrapper scripts to make Monte Carlo simulation
of Al2O3 on Al spectra easier.
Place this file in DTSA_ROOT/lib/dtsa2/    call with
import dtsa2.mc3AnodizedAl as mcAnodeAl"""


def simSonoraSpc(tAl2O3, e0, det, wkDst=5, lt=100, pc=1, nTraj=1000, xtraParams={}):
    """simSonoraSpc(tAl2O3, e0, det, wkDst=5, lt=100, pc=1, nTraj=1000, xtraParams={})
    Simulate a spectrum from tAl2O3 nm of Al2O3 on Al recoreed at
    e0 kV using the DTSA detector det and a wkDst mm working distance for
    lt sec with a probe current of pc nA. Compute nTraj trajectories."""
    tAl = 100 # um
    al2o3 = dt2.material("Al2O3",density=3.95)
    al  = dt2.material("Al", density=2.70)
    lAl2O3 = [al2o3, tAl2O3*1.0e-9]
    lAl = [al, tAl*1.0e-6]
    lay = [lAl2O3, lAl]
    sNam = "%g-nm-Al2O3-on-Al-%g-kV" % (tAl2O3, e0)
    spc = dt2.wrap(mc3.multiFilm(lay, det, e0, True, nTraj, lt*pc, True, True, xtraParams))
    props=spc.getProperties()
    props.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName, sNam)
    props.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
    props.setNumericProperty(epq.SpectrumProperties.FaradayBegin, pc)
    props.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
    props.setNumericProperty(epq.SpectrumProperties.WorkingDistance, wkDst)
    return(spc)

# Let's test
gitDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/"
simDir = gitDir + relPrj + "/mc3AnodizedAl/"
jmg.ensureDir(simDir)

imgSize =    512         # pixel size for images
imgSzUm =      5.0     # image size in microns
vmrlEl  =    100         # number of el for VMRL

wd = simDir
os.chdir(wd)
pyrDir = wd + "/sim-multifilm-on-sub Results"

det    = findDetector("Oxford p4 05eV 2K")
print(det)

xrts=mc3.suggestTransitions("AlO")
print(xrts)
for tr in xrts:
    print(tr.getSiegbahnName())
for tr in xrts:
    print(tr.getIUPACName())
# set up the extra parameters
xtraParams={}
xtraParams.update(mc3.configureXRayAccumulators(xrts, charAccum=True, charFluorAccum=True, bremFluorAccum=True))
# note that the image size on the specimen is in meters...
xtraParams.update(mc3.configureEmissionImages(xrts, imgSzUm*1.0e-6, imgSize))
xtraParams.update(mc3.configurePhiRhoZ(imgSzUm*1.0e-6))
xtraParams.update(mc3.configureTrajectoryImage(imgSzUm*1.0e-6, imgSize))
xtraParams.update(mc3.configureVRML(nElectrons = vmrlEl))
xtraParams.update(mc3.configureOutput(simDir))

print(xtraParams)

spc = simSonoraSpc(500, 30.0, det, xtraParams=xtraParams)
display(spc)