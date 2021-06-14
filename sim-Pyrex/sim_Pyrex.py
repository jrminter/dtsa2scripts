# -*- coding: utf-8 -*-
#         1         2         3         4         5         6         7 |
# 23456789012345678901234567890123456789012345678901234567890123456789012
#
# sim-Pyrex.py
# jrm 2021-06-07  starting a run with 1,000 traj
#                 Elapse: 
import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3, FluorescenceXRayGeneration3, XRayTransport3", None)
import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQLibrary.Detector as epd
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.microanalysis.NISTMonte.Gen3 as nm3
import gov.nist.microanalysis.EPQTools as et
import dtsa2.mcSimulate3 as mc3
import dtsa2.jmGen as jmg
import dtsa2.jmMC3 as jm3
import java.util as jutil
import java.io as jio
import java.nio.charset as cs
import os
import shutil
import time
import datetime


det      = findDetector("Oxford p4 05eV 2K")
e0       =    20    # kV
nTraj    = 50000    # trajectories
lt       =   500    # sec
pc       =     5.0  # nA
dose     = pc * lt  # na-sec"
bSaveSpc = True

print e0

homDir = os.environ['HOME']

print(homDir)

std = defineStd(("B", "O",   "Na",  "Mg", "Al", "Si",  "Cl",  "Ca", "Fe",),
                (3.92, 53.85, 3.12, 0.03, 1.17, 37.70, 0.10, 0.070, 0.030,),"Pyrex", 2.50)
                
###
          
def sim_amc_coated_mat(mat, det, e0, nTraj, lt=100, pc=1.0, tc=20.0):
    """sim_amc_coated_mat(mat, det, e0, nTraj, lt=100, pc=1.0, tc=20.0)

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
    tc - float (20.0)
        C thickness in nm


    Returns
    -------
    sim - DTSA scriptable spectrum 
        The simulated spectrum
    
    Example
    -------
    import dtsa2.jmMC3 as jm3
    det = findDetector("Oxford p4 05eV 2K")
    si = material("Si", density=2.3296)
    a = jm3.simCarbonCoatedStd(mgo, det, 5.0, 100, 100, 1.0, 20.0)
    a.display()

    """
    dose = pc * lt  # na-sec"
    amc = material("C", density=1.35)
    
    amcThickComment = "amC Thickness = %g nm %g trajectories" % (tc, nTraj)
    layers = [ [amc, tc*1.0e-9],
               [mat, 50.0e-6]
             ]
    xrts = []

    trs = mc3.suggestTransitions(amc, e0)
    for tr in trs:
       xrts.append(tr)

    trs = mc3.suggestTransitions(mat, e0)
    for tr in trs:
        xrts.append(tr)

    xtraParams={}
    xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))

    sim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                        dose=dose, sf=True, bf=True, xtraParams=xtraParams)
    sName = "%g-nm-amC-on-%s-%g-kV" % (tc, mat, e0)
    sim.rename(sName)
    sim.getProperties().setTextProperty(epq.SpectrumProperties.SpectrumComment, amcThickComment)
    return sim

start = time.time()

projDir = homDir + "\\OneDrive\\Documents\\git\\dtsa2scripts\\sim-Pyrex\\"
print(projDir)
jmg.ensureDir(projDir)

spc = sim_amc_coated_mat(std, det, e0, nTraj, lt=100, pc=1.0, tc=20.0)
display(spc)


# rptDir = homDir + relPrj + "/sim-C-on-Si-7kV Results"
# csvFil = homDir + relPrj + "/dtsa2-C-on-Si-%g-kV-kratios-%g-traj.csv" % (e0, nTraj)

###