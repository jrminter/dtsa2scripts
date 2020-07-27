
import sys
import os
import time
import shutil

# import dtsa2 as dt2
import dtsa2.mcSimulate3 as mc3

gitDir  = os.environ['GIT_HOME']
relPrj  = "/dtsa2Scripts/utility"
prjDir  = gitDir + relPrj
rptDir  = prjDir + '/testMultiBlock Results/'

det      = findDetector("Oxford p4 05eV 4K")
e0       =    15.0   # kV
nTraj    =  2000    # trajectories
lt       =  1000    # sec
pc       =     5.0   # nA
imgSize  =   512     # pixel size for images
imgSzUm  =    1.0    # image size in microns
vmrlEl   =    100     # number of el for VMRL
przDepUm =     5.0   # phirhoz depth in microns

start = time.time()

dose = pc*lt


nbsNiFe = mixture({"Si" : 0.0032,
                   "Cr" : 0.0006,
                   "Mn" : 0.0030,
                   "Fe" : 0.5100,
                   "Co" : 0.0002,
                   "Ni" : 0.4820,
                   "Cu" : 0.0004,
                   "Mo" : 0.0001},
                   8.2,
                   'nbsNiFe')
Epon828 =  mixture({"H"  : 0.0706,
                    "C"  : 0.7440,
                    "O"  : 0.1855,
                    "Cl" : 0.0030},
                     1.16,
                    'Epon828')

det = findDetector("Probe")

# use a limited set for emission images
xrtsEI = [transition("C K-L3"),
          transition("Fe K-L3"), transition("Fe L3-M5"),
          transition("Ni K-L3"), transition("Ni L3-M5")]

xrts = []

trs = mc3.suggestTransitions(nbsNiFe, e0)
for tr in trs:
    xrts.append(tr)

trs = mc3.suggestTransitions(Epon828, e0)
for tr in trs:
    xrts.append(tr)

blocks = [[[10.0e-6, 10.0e-6,10.0e-6], [-5.0e-6, 0, 1e-6], nbsNiFe],
          [[10.0e-6, 10.0e-6,10.0e-6], [5.0e-6, 0, 1e-6], Epon828]]

xtraParams={}
xtraParams.update(mc3.configurePhiRhoZ(przDepUm*1.0e-6))
xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))
# note that the image size on the specimen is in meters...
xtraParams.update(mc3.configureEmissionImages(xrtsEI, imgSzUm*1.0e-6, imgSize))
xtraParams.update(mc3.configureTrajectoryImage(imgSzUm*1.0e-6, imgSize))
xtraParams.update(mc3.configureVRML(nElectrons = vmrlEl))
# xtraParams.update(mc3.configureOutput(simDir))

print(xtraParams)


spc = mc3.multiblock(blocks, det, e0, True, nTraj, dose, True, True, xtraParams)

display(spc)



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