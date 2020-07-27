# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

measureCuL.py

A reproducible script to measure the bkg-subtracted CuL intensity
and uncertainty. This can be a proxy for probe current.

  Date      Who   Ver    Comment
----------  ---  -----   ------------------------------------------
2020-06-18  JRM  0.0.1   Initial test

"""

__revision__ = "$Id: measure_CuL.py John R. Minter $"
__version__ = "0.0.1"

import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3,BremsstrahlungXRayGeneration3,FluorescenceXRayGeneration3, XRayTransport3", None)

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


DataManager.clearSpectrumList()

start = time.time()

cu_spc_path = "D:/edsSpc/jrmOxfordSpc/03_XCS8-2017-12-14/Cu-2017-12-14-5eV-ch-2K.msa"
det = findDetector("Oxford p4 05eV 2K")
cu_spc = readSpectrum(cu_spc_path, 0, det)
cu_spc.display()

pi = cu_spc.peakIntegral(640, 1190) # now convert to cps

print("mean: %.2f cps, uncertainty: %.2f cps" % (pi.doubleValue()/60.0, pi.uncertainty()/60.0))

print(dir(pi))

# print(pi.mean(0.01))
# print(pi.uncertainty())


end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg





