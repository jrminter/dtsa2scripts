# -*- coding: utf-8 -*-
"""
DTSA-II Script - J. R. Minter - 2016-10-12

massFractionsTheEasyWay.py

  Date      Who  Comment
----------  ---  -----------------------------------------------
2016-10-12  JRM  Mass fractions the easy way...

Elapse: 0:00:00.0  ROCPW7ZC5C42

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
rptDir  = prjDir + '/massFractionsTheEasyWay Results/'

azo = material("Al2Zn98O100", density=5.61)

wfAl = round(azo.weightFractionU(epq.Element.Al, True).doubleValue(), 5)
wfZn = round(azo.weightFractionU(epq.Element.Zn, True).doubleValue(), 5)
wfO  = round(azo.weightFractionU(epq.Element.O, True).doubleValue(), 5)

AZO = {"Al" : wfAl, "Zn": wfZn, "O" : wfO}

print(AZO)

es = azo.getElementSet()

# print(dir(es))

# clean up cruft
shutil.rmtree(rptDir)
print "Done!"

