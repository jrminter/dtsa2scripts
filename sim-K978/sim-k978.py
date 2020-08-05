# -*- coding: utf-8 -*-
"""
sim-k978.py

"""

import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3,FluorescenceXRayGeneration3, XRayTransport3", None)


import os

import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQLibrary.Detector as epd
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.microanalysis.NISTMonte.Gen3 as nm3
import gov.nist.microanalysis.EPQTools as et
import gov.nist.microanalysis.dtsa2 as dt2

import dtsa2.jmGen as jmg
import dtsa2.mcSimulate3 as mc3
import dtsa2.hyperTools as ht

import shutil
import time

import java.util as jutil
import java.io as jio
import java.nio.charset as cs
import string


start = time.time()

outPath = "C:/Users/johnr/Documents/work/dtsa2-test-scripts/2020-08-07-sim-k978-glass"
det     = findDetector("SiLi")
e0      = 15      # kV
nTraj   = 1000  # trajectories
lt      = 60      # sec
pc      = 0.25    # nA
tNmC    = 20
dose = lt*pc

spcDir = outPath + "/sim-k978"
jmg.ensureDir(spcDir)
os.chdir(spcDir)


l_comps = [ epq.Element.Mg, epq.Element.Al, epq.Element.Si,
            epq.Element.K,  epq.Element.Ca, epq.Element.Ti,
            epq.Element.Rb, epq.Element.Sr, epq.Element.Cs,
            epq.Element.La, epq.Element.Eu, epq.Element.Gd,
            epq.Element.O]
l_mfs   = [ 0.013569,        0.089228,        0.21072,
            0.018678,        0.067467,        0.006714,
            0.020574,        0.100541,        0.002122,
            0.019185,        0.019432,        0.019520,
            0.393897]
k978 = jmg.create_material_from_mf(l_comps, l_mfs, 3.929, "K-978")

c = material("C", density=2.1)
layers = [ [c,    2.0e-9],
           [k978, 50.0e-6]
         ]

k978Sim = mc3.multiFilm(layers, det, e0, withPoisson=False, nTraj=nTraj,
                        dose=dose, sf=True, bf=True, xtraParams={})
sName = "%g-nm-C-on-K978" % tNmC
k978Sim.rename(sName)
k978Sim.setAsStandard(k978)
k978Sim.display()
fi =  spcDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (nTraj)
k978Sim.save(fi)

# clean up cruft
# shutil.rmtree(rptDir)
print "Done!"

end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg



