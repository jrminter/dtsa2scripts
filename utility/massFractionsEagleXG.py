# -*- coding: utf-8 -*-
"""
DTSA-II Script - J. R. Minter - 2018-07-30

massFractionsEagleXG.py

  Date      Who  Comment
----------  ---  -----------------------------------------------
2018-07-30  JRM  Mass fractions the easy way for EagleXG rev sort
2018-10-02  JRM  Change name to remove spaces. This made it easier
                 to add to the database.

Done!
Elapse: 0:00:00.6 jrmFastMac

 Z   Sym     WF
 8    O:   0.51877
14   Si:   0.30139
13   Al:   0.09008
20   Ca:   0.03876
 5    B:   0.03266
12   Mg:   0.00773
38   Sr:   0.00696
50   Sn:   0.00120
56   Ba:   0.00109
51   Sb:   0.00039
26   Fe:   0.00036
33   As:   0.00024
22   Ti:   0.00015

density = 2.36
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
rptDir  = prjDir + '/massFractionsEagleXG Results/'

eagleXG = mixture({"SiO2"  : 0.6447825,
                  "Al2O3" : 0.1702057,
                  "B2O3"  : 0.1051482,
                  "CaO"   : 0.0542376,
                  "MgO"   : 0.0128153,
                  "SrO"   : 0.0082368,
                  "SnO2"  : 0.0015215,
                  "BaO"   : 0.0012188,
                  "Fe2O3" : 0.0005078,
                  "Sb2O3" : 0.0004635,
                  "As2O3" : 0.0003145,
                  "ZrO2"  : 0.0002938,
                  "TiO2"  : 0.0002540
                 },
                 density=2.36,
                 name="eagleXG")

wfO  = round(eagleXG.weightFractionU(epq.Element.O,  True).doubleValue(), 5)
wfSi = round(eagleXG.weightFractionU(epq.Element.Si, True).doubleValue(), 5)
wfAl = round(eagleXG.weightFractionU(epq.Element.Al, True).doubleValue(), 5)
wfB  = round(eagleXG.weightFractionU(epq.Element.B,  True).doubleValue(), 5)
wfCa = round(eagleXG.weightFractionU(epq.Element.Ca, True).doubleValue(), 5)
wfMg = round(eagleXG.weightFractionU(epq.Element.Mg, True).doubleValue(), 5)
wfSr = round(eagleXG.weightFractionU(epq.Element.Sr, True).doubleValue(), 5)
wfSn = round(eagleXG.weightFractionU(epq.Element.Sn, True).doubleValue(), 5)
wfBa = round(eagleXG.weightFractionU(epq.Element.Ba, True).doubleValue(), 5)
wfFe = round(eagleXG.weightFractionU(epq.Element.Fe, True).doubleValue(), 5)
wfSb = round(eagleXG.weightFractionU(epq.Element.Sb, True).doubleValue(), 5)
wfAs = round(eagleXG.weightFractionU(epq.Element.As, True).doubleValue(), 5)
wfZr = round(eagleXG.weightFractionU(epq.Element.Zr, True).doubleValue(), 5)
wfTi = round(eagleXG.weightFractionU(epq.Element.Ti, True).doubleValue(), 5)


exg = {  "O" : wfO, 
        "Si" : wfSi, 
        "Al" : wfAl,
         "B" : wfB,
        "Ca" : wfCa,
        "Mg" : wfMg,
        "Sr" : wfSr,
        "Sn" : wfSn,
        "Ba" : wfBa,
        "Fe" : wfFe,
        "Sb" : wfSb,
        "As" : wfAs,
        "Ti" : wfTi
     }

# a1_sorted_keys = sorted(a1, key=a1.get, reverse=True)
# for r in a1_sorted_keys:
#    print r, a1[r]

for key, value in sorted(exg.iteritems(), key=lambda (k,v): (v,k), reverse=True):
    print("%s: %.5f" % (key, exg[key]))

# print(exg)

# es = exg.getElementSet()

# print(es)

# clean up cruft
shutil.rmtree(rptDir)
print "Done!"

