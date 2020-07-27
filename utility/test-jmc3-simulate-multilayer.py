"""
test-jmc3-simulate-multilayer.py

Extract emitted and detected spectra

Elapse: 0:00:08.3 for 100 traj
Elapse: 0:12:53.0 for 10000 traj for 250/100 nm
Elapse: 0:15:30.7 for 10000 traj for 500/250 nm

"""


import os, shutil, time
import dtsa2.jmcSimulate3 as jmc3
def ensureDir(d):
    """ensureDir(d)
    Check if the directory, d, exists, and if not create it."""
    if not os.path.exists(d):
        os.makedirs(d)

start = time.time()

gitDir = os.environ['GIT_HOME']
gitDir = gitDir.replace('\\','/')
wrkDir = gitDir + "/dtsa2Scripts/utility"
outDir = gitDir + "/dtsa2Scripts/utility/output/"
rptDir = wrkDir + '/test-jmc3-simulate-multilayer Results/'

DataManager.clearSpectrumList()

tNmIr   =   500.0  # nm of Ir top layer
tNmAg   =   250.0  # nm of Ag 2nd layer
nTraj   =  10000   # num Traj to run per pt 10000 for a long run
charF   =  True    # include characteristic fluorescence
bremF   =  True    # include continuum fluorescence 
pc      =     2.5  # nA
lt      =   100.0  # sec
e0      =    30.0  # keV

det = findDetector("Probe")
mat = material("Fe2O3",5.0)
dose    = pc * lt  # nA sec
titDet  = "Detected-%gnm-Ir-on-%g-nm-Ag-on-SiO2-%gkV" % (tNmIr, tNmAg, e0)
titEmi  = "Emitted-%gnm-Ir-on-%g-nm-Ag-on-SiO2-%gkV" % (tNmIr, tNmAg, e0)
# can use
# det  = findDetector("Si(Li)")
det  = findDetector("Oxford p4 05eV 4K")
print(det)

if 'defaultXtraParams' not in globals():
   defaultXtraParams = {}
if 'defaultBremFluor' not in globals():
   defaultBremFluor = False
if 'defaultCharFluor' not in globals():
   defaultCharFluor = False
if 'defaultNumTraj' not in globals():
   defaultNumTraj = 1000
if 'defaultDose' not in globals():
   defaultDose = 120.0

# start clean
DataManager.clearSpectrumList()
sio2 = material("SiO2", density= 2.65)
ir   = material("Ir",   density=22.56)
ag   = material("Ag",   density=11.9)

# define the desired transitions
xrts = jmc3.suggestTransitions("SiOIrAg")
xtraParams={}
xtraParams.update(jmc3.configureXRayAccumulators(xrts, charAccum=charF,
                  charFluorAccum=charF, bremFluorAccum=bremF))
print(xtraParams)

layers = [ [ir, tNmIr*1.0e-9],
           [ag, tNmAg*1.0e-9],
           [sio2, 50.0e-6]
         ]

lSpec = jmc3.multiFilm(layers, det, e0, True, nTraj,
                      dose, True, True, xtraParams)

# lSpec = jmc3.simulate(mat, det, e0, 120.0, True, 100, True, True, {})
lSpec[0].rename(titDet)
lSpec[0].display()
lSpec[1].rename(titEmi)
lSpec[1].display()

fi =  outDir + "/"
fi += titDet
fi += "-%g-Traj.msa" % (nTraj)
lSpec[0].save(fi)

fi =  outDir + "/"
fi += titEmi
fi += "-%g-Traj.msa" % (nTraj)
lSpec[1].save(fi)

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
