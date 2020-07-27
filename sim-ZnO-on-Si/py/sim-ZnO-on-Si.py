# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

sim-ZnO-on-Si.py

A script to simulate a thin film of 

  Date      Who   Ver    Comment
----------  ---  -----   ------------------------------------------
2019-12-19  JRM  0.0.1   Initial test. Required 38.168 min
                         

1> tabulate(selected())
Name  O        Al      Si       Ca       Ti     
K309  37.9945  7.8947  18.8471  10.7828  0.0527
      Fe       Zn      Ba       Total
      10.5587  0.8665  13.9387  100.9356

For some reason the residual looks like Bremsstrahlung...

"""

__revision__ = "$Id: sim-K309-quant.py John R. Minter $"
__version__ = "0.0.1"
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

def sim_zno_on_si(tzno, det, e0, nTraj, lt=100, pc=1.0):
    """sim_zno_on_si(tzno, det, e0, nTraj, lt=100, pc=1.0)

    Use mc3 multilayer simulation to simulate an am-C-ctd specimen

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
    lt - integer (100)
        The live time (sec)
    pc - float (1.0)
        The probe current in nA
    tzno - float (20.0)
        ZnO thickness in nm


    Returns
    -------
    sim - DTSA scriptable spectrum 
        The simulated spectrum
    
    Example
    -------
    import dtsa2.jmMC3 as jm3
    det = findDetector("Oxford p4 05eV 2K")
    si = material("Si", density=2.3296)
    a = jm3.simCarbonCoatedStd(mgo, det, 20.0, 100, 100, 1.0, 20.0)
    a.display()

    """
    dose = pc * lt  # na-sec"
    zno = material("ZnO", density=5.61)
    si = material("Si", density=2.329)

    znoThickComment = "ZnO Thickness = %g nm %g trajectories" % (tzno, nTraj)
    print(znoThickComment)
    layers = [ [zno, tzno*1.0e-9],
               [si, 50.0e-6]
             ]
    xrts = []

    trs = mc3.suggestTransitions(zno, e0)
    for tr in trs:
       xrts.append(tr)

    trs = mc3.suggestTransitions(si, e0)
    for tr in trs:
        xrts.append(tr)

    xtraParams={}
    xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))

    sim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                        dose=dose, sf=True, bf=True, xtraParams=xtraParams)
    sName = "%g-nm-ZnO-on-Si-%g-kV" % (tzno, e0)
    sim.rename(sName)
    sim.getProperties().setTextProperty(epq.SpectrumProperties.SpectrumComment, znoThickComment)
    return sim


nTraj   =  2000    # num Traj to run per pt 10000 for a long run
charF   =    True  # include characteristic fluorescence
bremF   =    True  # include continuum fluorescence 
pc      =    2.5   # nA
lt      =  100.0   # sec
e0      =    3.5   # keV
imgSize =  512     # pixel size for images
imgSzUm =   5.0    # image size in microns
vmrlEl  =   40     # number of el for VMRL
dose    = pc * lt  # nA sec

lNmZnOsim = [  1.0, 5.0, 10.0, 25.0, 50.0, 100.0]



gitDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/sim-ZnO-on-Si/"
simDir = gitDir + relPrj
jmg.ensureDir(simDir)

wd = gitDir + relPrj + "py"
os.chdir(wd)
pyrDir = wd + "/sim-multifilm-on-sub Results"

det  = findDetector("Oxford p4 05eV 2K")
print(det)


# start clean
DataManager.clearSpectrumList()
start = time.time()

zno = material("ZnO", density=5.61)
si  = material("Si", density=2.329)

znoSpc = mc3.simulate(zno, det, e0, dose, True, nTraj, True, True, {})
name = "ZnO at %g kV" % (e0)
znoSpc.rename(name)
znoSpc.setAsStandard(zno)
display(znoSpc)

siSpc = mc3.simulate(si, det, e0, dose, True, nTraj, True, True, {})
name = "Si at %g kV" % (e0)
siSpc.rename(name)
siSpc.setAsStandard(si)
display(siSpc)

i = 0
nSpec = len(lNmZnOsim)
for tzno in lNmZnOsim:
    i = i + 1
    startCycle = time.time()
    msg = "Simulating %g nm ZnO on Si at %g kV %g traj" % (tzno, e0, nTraj)
    print(msg)
    spc = sim_zno_on_si(tzno, det, e0, nTraj, lt=100, pc=1.0)
    display(spc)


end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg

