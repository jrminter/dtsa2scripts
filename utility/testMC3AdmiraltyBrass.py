"""
testMC3AdmiraltyBrass.py

Do a MC3 simulatiion of Admiralty Brass at 15 kV

"""

import os, shutil, time
import dtsa2.dt2Conv as dt2c
import dtsa2.mcSimulate3 as mc3

start = time.time()


"""
From:
text and here:
https://www.thoughtco.com/common-brass-alloys-and-their-uses-603706

densities
https://www.angstromsciences.com/density-elements-chart

The values in the dictionary are mass fractions of the component
oxides.
"""

name = 'Admiralty Brass'
rho = 0.69*8.86 + 0.30*7.13 + 0.01*7.31
mix = mixture({ "Cu"    : 0.69,
                "Zn"    : 0.30,
                "Sn"    : 0.01
               },
               density=rho,
               name=name)

det  = findDetector("Oxford p4 05eV 4K")
nTraj    = 200     # trajectories
lt       =   100     # sec
pc       =     1.0   # nA
e0       =    15.0   # kV
imgSize  =   512     # pixel size for images
imgSzUm  =    1.0    # image size in microns
vmrlEl   =    40     # number of el for VMRL
przDepUm =     3.0   # phirhoz depth in microns

dose = pc * lt  # na-sec"

# should not need to change below here

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
outDir = homDir + "/Documents/git/dtsa2Scripts/html"
rptDir = wrkDir + '/testMC3AdmiraltyBrass Results/'
simDir = homDir + "/Documents/git/dtsa2Scripts/sim-Admiralty-Brass"
dt2c.ensureDir(simDir)

DataManager.clearSpectrumList()

xrts = []
trs = mc3.suggestTransitions(mix, e0)
for tr in trs:
	xrts.append(tr)

# print(xrts)

xtraParams={}
xtraParams.update(mc3.configurePhiRhoZ(przDepUm*1.0e-6))
xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))
# note that the image size on the specimen is in meters...
xtraParams.update(mc3.configureEmissionImages(xrts, imgSzUm*1.0e-6,
                                              imgSize))
xtraParams.update(mc3.configureTrajectoryImage(imgSzUm*1.0e-6, imgSize))
xtraParams.update(mc3.configureVRML(nElectrons = vmrlEl))
xtraParams.update(mc3.configureOutput(simDir))

print(xtraParams)

sim = mc3.simulate(mix, det, e0, dose, True,
                   nTraj, True, True, xtraParams)
fmtS = "Admiralty-Brass-at-%g-kV"
sName = fmtS % (e0)
sim.rename(sName)
sim.display()
fi =  simDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (nTraj)
sim.save(fi)



# clean up cruft
shutil.rmtree(rptDir)
print "Done!"

end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg

