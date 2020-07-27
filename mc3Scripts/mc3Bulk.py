# mc3Bulk.py

# simple bulk simulation with specific transitons
# Date        Ver  Who  Notes
# 2015-10-01 0.90  JRM  Initial example. Verified with Iona v.2015-10-01

sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3, FluorescenceXRayGeneration3, XRayTransport3", None)
import dtsa2.mcSimulate3 as mc3
import os
import shutil
import dtsa2.jmGen as jmg


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
pyrDir = prjDir + "/mc3Bulk Results"


xrts = [transition("Fe K-L3"), transition("Fe K-M3"), transition("Fe L3-M5"), transition("O K-L3")]

xtraParams={}
xtraParams.update(mc3.configureXRayAccumulators(xrts, charAccum=True, charFluorAccum=True, bremFluorAccum=True))
# note that the image size on the specimen is in meters...
xtraParams.update(mc3.configureEmissionImages(xrts, imgSzUm*1.0e-6, imgSizePx))
xtraParams.update(mc3.configurePhiRhoZ(imgSzUm*1.0e-6))
xtraParams.update(mc3.configureTrajectoryImage(imgSzUm*1.0e-6, imgSizePx))
xtraParams.update(mc3.configureVRML(nElectrons = vmrlEl))
xtraParams.update(mc3.configureOutput(simDir))

spc = mc3.simulate(material("Fe2O3",5.0), det, e0, dose = pc*lt, nTraj=nTraj, sf=True, bf=True, xtraParams = xtraParams)

display(spc)

shutil.rmtree(pyrDir)
print "Done!"
