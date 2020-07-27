# testNewSimBulkStandard.py
#   Date      Who  Comment
# ----------  ---  -----------------------------------------------
# 2020-07-20  JRM  Try to have a more complete simBulkStd
#                  that can be inserted in a script
#                 

import dtsa2 as dtsa2
import dtsa2.mcSimulate3 as mc3
import time




def fullSimBulkStd(mat, det, e0, nTraj, outPath, dim=5.0e-6, lt=100, pc=1.0, emiSize=512, ctd=False):
    """
    fullSimBulkStd(mat, det, e0, nTraj, outPath, dim=5.0e-6, lt=100, pc=1.0, emiSize=512, ctd=False)

    Use mc3 simulation to simulate an uncoated standard specimen

    Parameters
    ----------
    mat - a dtsa material.
        Note the material must have an associated density. It should have a useful name.
    det - a dtsa detector
        Here is where we get the detector properties and calibration
    e0 - float
        The accelerating voltage in kV
    nTraj - integer
        The number of trajectories to run
    outPath - string
        The path to the directory for output
    dim - float (5.0e-6)
        The size of the emission images
    lt - integer (100)
        The live time (sec)
    pc - float (1.0)
        The probe current in nA
    emiSize - int (default 512)
        The width and depth of the emission images.
    ctd - Boolean (False) - is C coated


    Returns
    -------
    sim - DTSA scriptable spectrum 
        The simulated standard spectrum
    
    Example
    -------
    import dtsa2 as dtsa2
    import dtsa2.mcSimulate3 as mc3
    det = findDetector("Oxford p4 05eV 2K")
    cu = material("Cu", density=8.92)
    a = fullSimBulkStd(cu, det, 20.0, 100, 100, 1.0)
    a.display()

    """
    dose = pc * lt  # na-sec"
    xrts = []

    trs = mc3.suggestTransitions(mat, e0)
    for tr in trs:
        xrts.append(tr)
    mc3.configureEmissionImages(xrts, dim, emiSize)

    xtraParams={}
    xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))
    # note that the image size on the specimen is in meters...
    xtraParams.update(mc3.configureEmissionImages(xrts, 5.0e-6, 512))
    xtraParams.update(mc3.configurePhiRhoZ(5.0e-6))
    xtraParams.update(mc3.configureTrajectoryImage(5.0e-6, 512))
    xtraParams.update(mc3.configureVRML(nElectrons=100))
    xtraParams.update(mc3.configureOutput(outPath))
    mc3.configureOutput(outPath)
    print("Output sent to %s") % (outPath)

    dose = lt*pc

    sim = mc3.simulate(mat, det, e0, dose, withPoisson=True, nTraj=nTraj,
                       sf=True, bf=True, xtraParams=xtraParams)

    sName = "%s-%g-kV" % (mat, e0)
    sim.rename(sName)
    sim.setAsStandard(mat)
    sim.display()
    fi =  outPath + "/"
    fi += sName
    fi += "-%g-Traj.msa" % (nTraj)
    print(fi)
    sim.save(fi)
    return(sim)



det   = findDetector("Oxford p4 05eV 4K")
nTraj = 10000  # traj
dim   = 5.0e-6 # img dim microns
lt    = 100.0  # sec
pc    =  1.0   # nA
e0    = 20.0   # keV


cu = material("Cu", density=8.92)
ni = material("Ni", density=8.90)
fe = material("Fe", density=7.86)
ag = material("Ag", density=10.50)

DataManager.clearSpectrumList()

outPath = "C:/Users/johnr/Documents/work/test-dtsa-out"

start = time.time()

#        fullSimBulkStd(mat, det, e0,  nTraj, outPath, dim=5.0e-6, lt=100, pc=1.0, ctd=False)
specCu = fullSimBulkStd(cu,  det, e0,  nTraj, outPath, dim, lt, pc, False)

end = time.time()

delta = (end-start)/60
msg = "Cu simulation required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg

start = time.time()

specNi = fullSimBulkStd(ni,  det, e0,  nTraj, outPath, dim, lt, pc, False)

end = time.time()
delta = (end-start)/60
msg = "Ni simulation required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg


start = time.time()

specFe = fullSimBulkStd(fe,  det, e0,  nTraj, outPath, dim, lt, pc, False)

end = time.time()
delta = (end-start)/60
msg = "Fe simulation required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg

start = time.time()

specFe = fullSimBulkStd(fe,  det, e0,  nTraj, outPath, dim, lt, pc, False)

end = time.time()
delta = (end-start)/60
msg = "Fe simulation required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg

specAg = fullSimBulkStd(ag,  det, e0,  nTraj, outPath, dim, lt, pc, False)

end = time.time()
delta = (end-start)/60
msg = "Ag simulation required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg




# a = outPath + "/Cud-%g-nm-SiO2-on-C-%g-kV-%g-Traj.msa" % (tNmAuPd, tNmSiO2, e0, nTraj)
# a.save(fi)
