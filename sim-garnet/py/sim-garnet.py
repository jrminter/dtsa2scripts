# -*- coding: utf-8 -*-
# 
#         1         2         3         4         5         6         7 |
# 23456789012345678901234567890123456789012345678901234567890123456789012
#
# sim-garnet.py
# jrm 2019-12-20 - use mc3 to simulate garnet


import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3, FluorescenceXRayGeneration3, XRayTransport3", None)
import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQLibrary.Detector as epd
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.microanalysis.NISTMonte.Gen3 as nm3
import gov.nist.microanalysis.EPQTools as et
import dtsa2.mcSimulate3 as mc3
import dtsa2.jmGen as jmg
import java.util as jutil
import java.io as jio
import java.nio.charset as cs
import os
import shutil

def ensureDir(d):
  """ensureDir(d)
  Check if the directory, d, exists, and if not create it."""
  if not os.path.exists(d):
    os.makedirs(d)

gitDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/sim-garnet"
datDir = gitDir + relPrj + "/msa"
csvDir = gitDir + relPrj + "/csv"
ensureDir(csvDir)
wd = gitDir + relPrj + "/py"
os.chdir(wd)
pyrDir = wd + "/sim-garnet"
ensureDir(datDir)
det  = findDetector("Oxford p4 05eV 4K")
print(det)

nTraj    =  1000   # num Traj to run per pt 10000 for a long run
charF    =  True   # include characteristic fluorescence
bremF    =  True   # include continuum fluorescence 
pc       =     2.5  # nA
lt       =   100.0  # sec
e0       =    15.0  # keV

dose = pc * lt # nA sec
# start clean
DataManager.clearSpectrumList()


"""
Garnet 15 kV std
O  36.02
Mg  2.23
Al 11.55
Si 18.40
Ca  1.51
Mn 10.51
Fe 19.79

Calczaf

St 2524 Garnet
TakeOff = 40.0  KiloVolt = 15.0  Density =  5.000

A typical garnet composition. From a McCrone video on YouTube
Elemental Composition

Average Total Oxygen:         .000     Average Total Weight%:  100.010
Average Calculated Oxygen:    .000     Average Atomic Number:   15.300
Average Excess Oxygen:        .000     Average Atomic Weight:   24.943

ELEM:        O      Mg      Al      Si      Ca      Mn      Fe
XRAY:      ka      ka      ka      ka      ka      ka      ka 
ELWT:   36.020   2.230  11.550  18.400   1.510  10.510  19.790
KFAC:    .1934   .0132   .0797   .1364   .0145   .0915   .1743
ZCOR:   1.8621  1.6868  1.4499  1.3489  1.0384  1.1485  1.1353
AT% :   56.147   2.288  10.676  16.339    .940   4.771   8.838

"""


garnet = epq.Material(epq.Composition([epq.Element.O,
                                    epq.Element.Mg,
                                    epq.Element.Al,
                                    epq.Element.Si,
                                    epq.Element.Ca,
                                    epq.Element.Mn,
                                    epq.Element.Fe],
                                    [ 0.3602,
                                      0.0223,
                                      0.1155,
                                      0.1840,
                                      0.0151,
                                      0.1051,
                                      0.1979]), epq.ToSI.gPerCC(3.56))
garnet.setName("garnet")

trs = majorTransitions(garnet, e0, trs=(epq.XRayTransition.KA1, epq.XRayTransition.KB1, epq.XRayTransition.LA1, epq.XRayTransition.MA1), thresh=1.0)


garnetSpc = mc3.simulate(garnet, det, e0=e0, dose=dose, withPoisson=True, nTraj=nTraj, sf=charF, bf=bremF, xtraParams={})
sp=garnetSpc.getProperties()
sp.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName,"garent-std")
sp.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
sp.setNumericProperty(epq.SpectrumProperties.FaradayBegin, dose)
sp.setNumericProperty(epq.SpectrumProperties.FaradayEnd, dose)
sp.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
sp.setCompositionProperty(epq.SpectrumProperties.StandardComposition, epq.Composition(garnet))
display(garnetSpc)

"""

# clean up cruft
# shutil.rmtree(pyrDir)
print "Done!"




