# mc3CoatedSphereOnSubstrate.py

# Simulation of an Pt coated Al2O3 sphere on a C substrate
#
# Date        Ver  Who  Notes
# 2015-10-01 0.90  JRM  Initial example. Verified with Iona v.2015-10-01

sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3, FluorescenceXRayGeneration3, XRayTransport3", None)
import dtsa2.mcSimulate3 as mc3
import os
import shutil
import dtsa2.jmGen as jmg

radUm     =    0.5   # sphere radius in microns
ctgUm     =    0.010 # coating thickness
det       = findDetector("Probe") # DTSA-II default detector, use yours
e0        =   20     # keV
pc        =    2.5   # nA
lt        =  100.0   # sec
imgSzUm   =    5.0   # physical size of images in mucrons
imgSizePx =  512     # size of images in pixels
nTraj     = 1000     # number of trajectories
vmrlEl    =   40     # number of el for VMRL
dose      = pc * lt  # nA sec

gitHom = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/mc3Scripts"
prjDir = gitHom + relPrj
datDir = prjDir + "/dat"
simDir = datDir + "/sim"

jmg.ensureDir(datDir)
jmg.ensureDir(simDir)

os.chdir(prjDir)
pyrDir = prjDir + "/mc3CoatedSphereOnSubstrate Results"

al2o3 = epq.Material(epq.Composition([epq.Element.Al,epq.Element.O],[0.5293,0.4707]),epq.ToSI.gPerCC(3.95))
al2o3.setName("Al2O3")
c = epq.Material(epq.Composition([epq.Element.C],[1.0]),epq.ToSI.gPerCC(2.62))
c.setName("C")
pt = epq.Material(epq.Composition([epq.Element.Pt],[1.0]),epq.ToSI.gPerCC(21.4))
pt.setName("Pt")


xrts = [transition("Al K-L3"), transition("O K-L3"), transition("C K-L3"), transition("Pt M5-N7"), transition("Pt L3-M5")]

xtraParams={}
xtraParams.update(mc3.configureXRayAccumulators(xrts, charAccum=True, charFluorAccum=True, bremFluorAccum=True))
# note that the image size on the specimen is in meters...
xtraParams.update(mc3.configureEmissionImages(xrts, imgSzUm*1.0e-6, imgSizePx))
xtraParams.update(mc3.configurePhiRhoZ(imgSzUm*1.0e-6))
xtraParams.update(mc3.configureTrajectoryImage(imgSzUm*1.0e-6, imgSizePx))
xtraParams.update(mc3.configureVRML(nElectrons = vmrlEl))
xtraParams.update(mc3.configureOutput(simDir))

spc = mc3.coatedSphere(al2o3, radUm*1.0e-6, pt, ctgUm*1.0e-6, det, e0=e0, withPoisson=True, nTraj=nTraj, dose=dose, sf=True, bf=True, substrate=c, xtraParams=xtraParams)


display(spc)

shutil.rmtree(pyrDir)
print "Done!"
