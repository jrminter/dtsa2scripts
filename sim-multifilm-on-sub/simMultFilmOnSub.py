
import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3, FluorescenceXRayGeneration3, XRayTransport3", None)
import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQLibrary.Detector as epd
import gov.nist.microanalysis.Utility as epu
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.microanalysis.NISTMonte.Gen3 as nm3
import gov.nist.microanalysis.EPQTools as et
import dtsa2.mcSimulate3 as mc3
import dtsa2.jmGen as jmg
import dtsa2 as dt2
import java.util as jutil
import java.io as jio
import java.nio.charset as cs
import os
import glob
import shutil

def multiLayersOnSubstrate(layers, sub, det, e0=20.0, pc=1.0, lt=120.,
                           nTraj=100, withPoisson=True, sf=True,
                           bf=True, xtraParams={}):

    """multiLayersOnSubstrate(layers, sub, det, e0=20.0, pc=1.0, lt=100.,
                              nTraj=100, withPoisson=True, sf=True,
                              bf=True, xtraParams={})

    Monte Carlo simulate a spectrum from a multilayer thin film on a
    substrate.

    Parameters
    ----------

    layers: iterable list
        Layers is a iterable list of [material, thickness (meters)].
        Note that the materials must have associated densities.
    sub: material
        The material to use for the substrate
    det: detector
        The detector to use
    e0: float (20.0)
        The accelerating voltage in kV
    pc: float (1.0)
        The probe current in nA
    lt: float (100.0)
        the spectrum live time (sec)
    nTraj: an integer (100)
        Number of trajectories to compute
    withPoisson: boolean (True)
        Flag to include Poisson noise
    sf: boolean (True)
        Flag to include spectrum fluorescence
    bf: boolean (True)
        Flag to include bremsstrahlung fluorescence
    xtraParams: dictionary ({})
        A dictionary of extra parameters
    
    Returns
    -------
    spc: the simulated spectrum

    """
    dose = pc * lt
    # sl = "%s" % (",".join("%0.2f um of %s" % (1.0e6 * layer[1], layer[0]) for layer in layers)
    tmp = u"MC simulation of a multilayer film" #  [%s] on %s at %0.1f keV%s%s" %(sl , sub, e0, (" + CSF" if sf else ""), (" + BSF" if bf else ""))

    def buildFilm(monte, chamber, origin, buildParams):
        sr = chamber
        sub = buildParams["Sub"]
        for layer in buildParams["Layers"]:
             if layer[1] <= 0.0:
                raise "The layer thickness must be larger than zero."
             sr = monte.addSubRegion(sr, layer[0], nm.MultiPlaneShape.createSubstrate([0.0, 0.0, -1.0], origin))
             origin = epu.Math2.plus(origin, [0.0, 0.0, layer[1]])
        sr = monte.addSubRegion(sr, sub, nm.MultiPlaneShape.createSubstrate([0.0, 0.0, -1.0], origin))
    
    
    return mc3.base(det, e0, withPoisson, nTraj, dose, sf, bf, tmp, buildFilm, {"Layers": layers , "Sub": sub}, xtraParams)



widNm   =   25.0   # width of layers [nm]
nTraj   =  250     # num Traj to run per pt 10000 for a long run
charF   =    True  # include characteristic fluorescence
bremF   =    True  # include continuum fluorescence 
pc      =    2.5   # nA
lt      =  100.0   # sec
e0      =    7.0   # keV
imgSize =  512     # pixel size for images
imgSzUm =   5.0   # image size in microns
vmrlEl  =   40     # number of el for VMRL
dose    = pc * lt # nA sec


gitDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/"
simDir = gitDir + relPrj + "/"
jmg.ensureDir(simDir)

wd = gitDir + relPrj 
os.chdir(wd)
pyrDir = wd + "/sim-multifilm-on-sub Results"

det  = findDetector("Oxford p4 05eV 2K")
print(det)


# start clean
DataManager.clearSpectrumList()

al2o3 = epq.Material(epq.Composition([epq.Element.Al,epq.Element.O],[0.5293,0.4707]),epq.ToSI.gPerCC(3.95))
al = epq.Material(epq.Composition([epq.Element.Al],[1.0],"Al"), epq.ToSI.gPerCC(2.7))