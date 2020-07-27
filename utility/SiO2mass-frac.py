# -*- coding: utf-8 -*-
"""
DTSA-II Script - J. R. Minter - 2016-10-12

SiO2mass-frac.py

  Date      Who  Comment
----------  ---  -----------------------------------------------
2018-07-18  JRM  SiO2 mass fractions the easy way...

Elapse: 0:00:00.0  ROCPW7ZC5C42

Ir 1.0 22.56 g/cm   Z = 77
Ag 1.0 10.50 g/cm3  Z = 47
SiO2    2.65 g/cm3
Z = 8             14
  {'O': 0.53256, 'Si': 0.46744}

"""
import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3,FluorescenceXRayGeneration3, XRayTransport3", None)

import os
import glob
import shutil
import time
import math
import csv

import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.microanalysis.NISTMonte.Gen3 as nm3
import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQLibrary.Detector as epd
import gov.nist.microanalysis.Utility as epu
import gov.nist.microanalysis.EPQTools as ept


import dtsa2 as dt2
import dtsa2.mcSimulate3 as mc3

gitDir  = os.environ['GIT_HOME']
relPrj  = "/dtsa2Scripts/utility"
prjDir  = gitDir + relPrj
rptDir  = prjDir + '/SiO2mass-frac Results/'

sio2 = material("SiO2", density=2.65)

wfSi = round(sio2.weightFractionU(epq.Element.Si, True).doubleValue(), 5)
wfO  = round(sio2.weightFractionU(epq.Element.O, True).doubleValue(), 5)

Silica = { "Si": wfSi, "O" : wfO}

print(Silica)

es = sio2.getElementSet()

# print(dir(es))

# clean up cruft
shutil.rmtree(rptDir)
print "Done!"

