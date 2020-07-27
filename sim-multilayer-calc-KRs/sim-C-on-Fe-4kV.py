# -*- coding: utf-8 -*-
#         1         2         3         4         5         6         7 |
# 23456789012345678901234567890123456789012345678901234567890123456789012
#
# sim-C-on-Fe-4kV.py
# jrm 2018-09-22  starting a run with 500 traj
#.                This script required 1.955 min
#                 Elapse: 0:01:58.6
# jrm 2018-09-22  starting a run with 50,000 traj
#.                This script required 188.683 min
#                 Elapse: 3:08:41.7
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


det = findDetector("Oxford p4 05eV 2K")
e0       =     4    # kV
nTraj    = 50000    # trajectories
lt       =   500    # sec
pc       =     5.0 # nA
dose     = pc * lt  # na-sec"
bSaveSpc = True

lNmCsim = [  1.0,   5.0,  10.0,  15.0,  20.0,  25.0,  30.0,  35.0,  40.0,
            45.0,  50.0,  60.0,  70.0,  80.0,  90.0, 100.0, 110.0, 120.0, 
           130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0]


homDir = os.environ['HOME']
relPrj = "/Documents/work/dtsa2/sim-C-on-Fe-4kV"
datDir = homDir + relPrj + "/msa-%g" % (nTraj)
jmg.ensureDir(datDir)
rptDir = homDir + relPrj + "/sim-C-on-Fe-4kV Results"
csvFil = homDir + relPrj + "/dtsa2-C-on-Fe-%g-kV-kratios-%g-traj.csv" % (e0, nTraj)

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
    a = jm3.simCarbonCoatedStd(mgo, det, 20.0, 100, 100, 1.0, 20.0)
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

DataManager.clearSpectrumList()

c   = material("C", density=2.266)
fe  = material("Fe", density=7.874)

# only get FeL at 7 kV
trs = [epq.XRayTransitionSet(epq.Element.C, epq.XRayTransitionSet.K_FAMILY),
       epq.XRayTransitionSet(epq.Element.Fe, epq.XRayTransitionSet.L_FAMILY)]


# First simulate the standards
startCycle = time.time()
c_std_spc = jm3.simBulkStd(c, det, e0, nTraj, lt=lt, pc=pc)
display(c_std_spc)
sName = "%s-std-%g-kV-%g-traj" % ("C", e0, nTraj)
fi =  datDir + "/"
fi += sName
fi += ".msa"
if bSaveSpc == True:
    c_std_spc.save(fi)
endCycle = time.time()
delta = (endCycle-startCycle)/60
tod = datetime.datetime.now().time()
msg = "C std required %.3f min %s" % (delta, tod)
print msg

startCycle = time.time()
fe_std_spc = jm3.simBulkStd(fe, det, e0, nTraj, lt=lt, pc=pc)
display(fe_std_spc)
sName = "%s-std-%g-kV-%g-traj" % ("Fe", e0, nTraj)
fi =  datDir + "/"
fi += sName
fi += ".msa"
if bSaveSpc == True:
    fe_std_spc.save(fi)
endCycle = time.time()
delta = (endCycle-startCycle)/60
tod = datetime.datetime.now().time()
msg = "Fe std required %.3f min %s" % (delta, tod)
print msg

cStd  = {"El":element("C"),  "Spc":c_std_spc}
feStd = {"El":element("Fe"), "Spc":fe_std_spc}

stds = [cStd, feStd]

# seed the zero
lNmC = [0.0]
lkC  = [0.0]
lkFe = [1.0]

i = 0
nSpec = len(lNmCsim)
for tc in lNmCsim:
    i = i + 1
    startCycle = time.time()
    msg = "Simulating %g nm C on Fe at %g kV %g traj" % (tc, e0, nTraj)
    print(msg)
    spc = sim_amc_coated_mat(fe, det, e0, nTraj, lt=lt, pc=pc, tc=tc)
    display(spc)
    sName = "%g-nm-C-on-Fe-%g-kV-%g-traj" % (tc, e0, nTraj)
    fi =  datDir + "/"
    fi += sName
    fi += ".msa"
    if bSaveSpc == True:
        spc.save(fi)
    endCycle = time.time()
    delta = (endCycle-startCycle)/60
    tod = datetime.datetime.now().time()
    msg = "C layer %g of %g required %.3f min %s" % (i, nSpec, delta, tod)
    print msg
    res = jmg.compKRs(spc, stds, trs, det, e0)
    print(res)
    lNmC.append(tc)
    lkC.append(round(res[0], 5))
    lkFe.append(round(res[1], 5))

nMeas = len(lkFe)
f = open(csvFil, 'w')
strLine = 't.nm, kC, kFe\n'
f.write(strLine)
for i in range(nMeas):
    strLine = "%g" % lNmC[i] + ","
    strLine = strLine + "%.5f" % lkC[i] + ","
    strLine = strLine + "%.5f" % lkFe[i] + "\n"
    f.write(strLine)  
f.close()


# clean up cruft
shutil.rmtree(rptDir)
print "Done!"

end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg

