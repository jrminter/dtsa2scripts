# simulateCoatedK411Ex.py
#   Date      Who  Comment
# ----------  ---  -----------------------------------------------
# 2017-06-17  JRM  Simulate Ritchie K411 example
#                  Elapse: 1:40:39.7 on crunch

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

rptDir = wrkDir + '/simulateCoatedK411Ex Results/'

spcDir = wrkDir + "/dat/spc"
simDir = wrkDir + "/dat/sim"
jmg.ensureDir(spcDir)
jmg.ensureDir(simDir)
jmg.ensureDir(spcDir)

print(simDir)

det      = findDetector("Oxford p4 05eV 2K")
e0       =    20     # kV
nTraj    = 20000       # trajectories
lt       =   100     # sec
pc       =     1.0   # nA
tNmC     =    20     # nm C
przDepUm =     2.0   # depth for phi-rho-z images in microns
imgSzUm  =     2.0   # size for emission images in microns
imgSize  =   512     # pixel size for images
vmrlEl   =    40     # number of el for VMRL
resln    =     1.0   # factor for phirhoz resln


dose = pc * lt  # na-sec"

DataManager.clearSpectrumList()

apatite = epq.Material(epq.Composition([epq.Element.O, epq.Element.Ca, epq.Element.P, epq.Element.F, epq.Element.Sr],
                                       [0.39368,       0.38936,        0.17863,       0.03700,       0.00133]),
                                       epq.ToSI.gPerCC( 3.15))
apatite.setName("apatite")

k411 = epq.Material(epq.Composition([epq.Element.O, epq.Element.Si, epq.Element.Fe, epq.Element.Ca, epq.Element.Mg],
                                    [0.435581,      0.25382,        0.11209,        0.11057,        0.08847]),
                                    epq.ToSI.gPerCC(2.6))
k411.setName("K411")

mgo = material("MgO", density=3.58)
c = material("C", density=2.1)
si = material("Si", density=2.329)
fe = material("Fe", density= 7.86)
caf2 = material("CaF2", density= 3.18)

# start with the 'unknown'

layers = [ [c, tNmC*1.0e-9],
           [k411, 50.0e-6]
         ]

# multiFilm(layers, det, e0=20.0, withPoisson=True,
#           nTraj=defaultNumTraj, dose=defaultDose,
#           sf=defaultCharFluor, bf=defaultBremFluor,
#            xtraParams=defaultXtraParams)

xrts = []

trs = mc3.suggestTransitions(c, e0)
for tr in trs:
    xrts.append(tr)

trs = mc3.suggestTransitions(k411, e0)
for tr in trs:
    xrts.append(tr)

xtraParams={}
xtraParams.update(mc3.configurePhiRhoZ(przDepUm*1.0e-6, resln))
xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))
# note that the image size on the specimen is in meters...
xtraParams.update(mc3.configureEmissionImages(xrts, imgSzUm*1.0e-6, imgSize))
xtraParams.update(mc3.configureTrajectoryImage(imgSzUm*1.0e-6, imgSize))
xtraParams.update(mc3.configureVRML(nElectrons = vmrlEl))
xtraParams.update(mc3.configureOutput(simDir))

print(xtraParams)



k411Sim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                        dose=dose, sf=True, bf=True, xtraParams=xtraParams)
sName = "%g-nm-C-on-K411-%g-kV" % (tNmC, e0)
k411Sim.rename(sName)
k411Sim.setAsStandard(k411)
k411Sim.display()

fi = spcDir + "/" + sName + "-%g-Traj.msa" % (nTraj)
k411Sim.save(fi)

# Now apatite std

layers = [ [c, tNmC*1.0e-9],
           [apatite, 50.0e-6]
         ]

xrts = []

trs = mc3.suggestTransitions(c, e0)
for tr in trs:
    xrts.append(tr)

trs = mc3.suggestTransitions(apatite, e0)
for tr in trs:
    xrts.append(tr)

xtraParams={}
xtraParams.update(mc3.configurePhiRhoZ(przDepUm*1.0e-6, resln))
xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))
# note that the image size on the specimen is in meters...
xtraParams.update(mc3.configureEmissionImages(xrts, imgSzUm*1.0e-6, imgSize))
xtraParams.update(mc3.configureTrajectoryImage(imgSzUm*1.0e-6, imgSize))
xtraParams.update(mc3.configureVRML(nElectrons = vmrlEl))
xtraParams.update(mc3.configureOutput(simDir))

print(xtraParams)

apatiteSim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                           dose=dose, sf=True, bf=True, xtraParams=xtraParams)
sName = "%g-nm-C-on-apatite-%g-kV" % (tNmC, e0)
apatiteSim.rename(sName)
apatiteSim.setAsStandard(apatite)
apatiteSim.display()

fi = spcDir + "/" + sName + "-%g-Traj.msa" % (nTraj)
apatiteSim.save(fi)

# Now MgO std

layers = [ [c, tNmC*1.0e-9],
           [mgo, 50.0e-6]
         ]

xrts = []

trs = mc3.suggestTransitions(c, e0)
for tr in trs:
    xrts.append(tr)

trs = mc3.suggestTransitions(mgo, e0)
for tr in trs:
    xrts.append(tr)

xtraParams={}
xtraParams.update(mc3.configurePhiRhoZ(przDepUm*1.0e-6, resln))
xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))
# note that the image size on the specimen is in meters...
xtraParams.update(mc3.configureEmissionImages(xrts, imgSzUm*1.0e-6, imgSize))
xtraParams.update(mc3.configureTrajectoryImage(imgSzUm*1.0e-6, imgSize))
xtraParams.update(mc3.configureVRML(nElectrons = vmrlEl))
xtraParams.update(mc3.configureOutput(simDir))

print(xtraParams)

mgoSim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                       dose=dose, sf=True, bf=True, xtraParams=xtraParams)
sName = "%g-nm-C-on-MgO-%g-kV" % (tNmC, e0)
mgoSim.rename(sName)
mgoSim.setAsStandard(mgo)
mgoSim.display()

fi = spcDir + "/" + sName + "-%g-Traj.msa" % (nTraj)
mgoSim.save(fi)

# Now Si std

layers = [ [c, tNmC*1.0e-9],
           [si, 50.0e-6]
         ]

xrts = []

trs = mc3.suggestTransitions(c, e0)
for tr in trs:
    xrts.append(tr)

trs = mc3.suggestTransitions(si, e0)
for tr in trs:
    xrts.append(tr)

xtraParams={}
xtraParams.update(mc3.configurePhiRhoZ(przDepUm*1.0e-6, resln))
xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))
# note that the image size on the specimen is in meters...
xtraParams.update(mc3.configureEmissionImages(xrts, imgSzUm*1.0e-6, imgSize))
xtraParams.update(mc3.configureTrajectoryImage(imgSzUm*1.0e-6, imgSize))
xtraParams.update(mc3.configureVRML(nElectrons = vmrlEl))
xtraParams.update(mc3.configureOutput(simDir))

print(xtraParams)

siSim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                      dose=dose, sf=True, bf=True, xtraParams=xtraParams)
sName = "%g-nm-C-on-Si-%g-kV" % (tNmC, e0)
siSim.rename(sName)
siSim.setAsStandard(si)
siSim.display()

fi = spcDir + "/" + sName + "-%g-Traj.msa" % (nTraj)
siSim.save(fi)

# Now Fe std

layers = [ [c, tNmC*1.0e-9],
           [fe, 50.0e-6]
         ]

xrts = []

trs = mc3.suggestTransitions(c, e0)
for tr in trs:
    xrts.append(tr)

trs = mc3.suggestTransitions(fe, e0)
for tr in trs:
    xrts.append(tr)

xtraParams={}
xtraParams.update(mc3.configurePhiRhoZ(przDepUm*1.0e-6, resln))
xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))
# note that the image size on the specimen is in meters...
xtraParams.update(mc3.configureEmissionImages(xrts, imgSzUm*1.0e-6, imgSize))
xtraParams.update(mc3.configureTrajectoryImage(imgSzUm*1.0e-6, imgSize))
xtraParams.update(mc3.configureVRML(nElectrons = vmrlEl))
xtraParams.update(mc3.configureOutput(simDir))

print(xtraParams)

feSim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                      dose=dose, sf=True, bf=True, xtraParams=xtraParams)
sName = "%g-nm-C-on-Fe-%g-kV" % (tNmC, e0)
feSim.rename(sName)
feSim.setAsStandard(fe)
feSim.display()

fi = spcDir + "/" + sName + "-%g-Traj.msa" % (nTraj)
feSim.save(fi)

# Now CaF2 ref

layers = [ [c, tNmC*1.0e-9],
           [caf2, 50.0e-6]
         ]

xrts = []

trs = mc3.suggestTransitions(c, e0)
for tr in trs:
    xrts.append(tr)

trs = mc3.suggestTransitions(caf2, e0)
for tr in trs:
    xrts.append(tr)

xtraParams={}
xtraParams.update(mc3.configurePhiRhoZ(przDepUm*1.0e-6, resln))
xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))
# note that the image size on the specimen is in meters...
xtraParams.update(mc3.configureEmissionImages(xrts, imgSzUm*1.0e-6, imgSize))
xtraParams.update(mc3.configureTrajectoryImage(imgSzUm*1.0e-6, imgSize))
xtraParams.update(mc3.configureVRML(nElectrons = vmrlEl))
xtraParams.update(mc3.configureOutput(simDir))

print(xtraParams)

caf2Sim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                        dose=dose, sf=True, bf=True, xtraParams=xtraParams)
sName = "%g-nm-C-on-CaF2-%g-kV" % (tNmC, e0)
caf2Sim.rename(sName)
caf2Sim.setAsStandard(caf2)
caf2Sim.display()

fi = spcDir + "/" + sName + "-%g-Traj.msa" % (nTraj)
caf2Sim.save(fi)


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


