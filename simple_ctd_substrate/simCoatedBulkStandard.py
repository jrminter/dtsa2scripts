import dtsa2 as dtsa2
import dtsa2.mcSimulate3 as mc3
import dtsa2.jmMC3 

DataManager.clearSpectrumList()


outPath = "C:/Users/johnr/Documents/work/dtsa2-test-scripts/simCoatedBulkStandard"
e0      = 15 # accelerating voltage
# det = findDetector("Si(Li)") # Replace with your detector's name
det     = findDetector("Oxford p4 05eV 4K") # Replace with your detector's name
nTraj   = 10000   # trajectories
lt      = 60      # sec
pc      = 0.25    # nA
tNmC    = 20

test = {}


ctg = material("C", 1.8)
mat = material("Cu", 8.96)

def fullSimBulkStd(mat, ctg, ctgThickNm, det, e0, nTraj, outPath,
                   dim=5.0e-6, lt=100, pc=1.0, emiSize=512, ctd=False):
    """
    fullSimBulkStd(mat, ctg, ctgThickNm, det, e0, nTraj, outPath,
                   dim=5.0e-6, lt=100, pc=1.0, emiSize=512, ctd=False)

    Use mc3 simulation to simulate an uncoated standard specimen

    Parameters
    ----------
    mat - a dtsa material.
        Note the material must have an associated density. It should
        have a useful name.

    ctg - a dtsa2 material for the coating

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
    det = findDetector("Oxford p4 05eV 4K")
    cu = material("Cu", density=8.92)

    outPath = "C:/Users/johnr/Documents/git/dtsa2Scripts/ben-buse/out"

    a = fullSimBulkStd(cu, det, 15.0, 100, outPath,
                       dim=5.0e-6, lt=100,
                       pc=1.0, emiSize=512, ctd=False)
    a.display()

    """
    strCtg = ctg.getName()
    strMat = mat.getName()
    dose = pc * lt  # na-sec"
    
    # specify the transitions to generate
    xrts = []

    trs = mc3.suggestTransitions(mat, e0)
    for tr in trs:
        xrts.append(tr)
    trs = mc3.suggestTransitions(ctg, e0)
    for tr in trs:
        xrts.append(tr)

    # At 20 kV the images are best at 2.0e-6
    xtraParams={}
    xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))
    # note that the image size on the specimen is in meters...
    xtraParams.update(mc3.configureEmissionImages(xrts, 2.0e-6, 512))
    xtraParams.update(mc3.configurePhiRhoZ(2.0e-6))
    xtraParams.update(mc3.configureTrajectoryImage(2.0e-6, 512))
    xtraParams.update(mc3.configureVRML(nElectrons=100))
    xtraParams.update(mc3.configureOutput(outPath))
    print("Output sent to %s") % (outPath)
    layers = [ [ctg, ctgThickNm*1.0e-9],
               [mat, 1.0e-3]]
    sim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                        dose=dose, sf=True, bf=True, xtraParams=xtraParams)
    sName = "%g-nm-%s-on-%s" % (tNmC, strCtg, strMat)
    sim.rename(sName)
    sim.setAsStandard(mat)
    sim.display()
    return(sim)



out = fullSimBulkStd(mat, ctg, tNmC, det, e0, nTraj, outPath,
                     dim=2.0e-6, lt=100, pc=1.0, emiSize=512, ctd=False)
out.display()

strCtg = ctg.getName()
strMat = mat.getName()

fi =  outPath + "/"
fi += "%g-nm-%s-on-%s" % (tNmC, strCtg, strMat)
fi += "-%g-Traj.msa" % (nTraj)
out.save(fi)