# anaDet-P4-5eV-4K.py
#
#   Date      Who  Comment
# ----------  ---  -----------------------------------------------
# 2016-06-17  JRM  Analyze a series of Cu spectra from the Sirion
#                  to characterize 'Oxford p4 05eV 2K'

import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3,FluorescenceXRayGeneration3, XRayTransport3", None)


import os
import glob

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

srcDir = os.environ['GIT_HOME']
prjDir = "/dtsa2Scripts/anaSirion/py/"
relDat = "/dtsa2Scripts/anaSirion/dat/"
relSpc = "/dtsa2Scripts/spc/Sirion/Oxford-P4-05eV-4K"
relPy  = "anaDet-P4-5eV-4K Results/"
rptDir = srcDir + prjDir + relPy
datDir = srcDir + relDat
spcDir = srcDir + relSpc

det = findDetector("Oxford p4 05eV 4K")
mIterations = 5
e0 = 15
cuStd = Database.findStandard("Cu standard")
bDisplayAll = False


csvFil = "%s/Oxford-P4-05eV-4K-%gkV-Cu.csv" % (datDir, e0)
print(csvFil)

query = spcDir + '*.msa'
print(query)

lDate = []
lCwMu = []
lCwUN = []
lZoMu = []
lZoUN = []
lRsMu = []
lRsUN = []
lPkMu = []
lPkUN = []

for name in glob.glob(spcDir + '/*.msa'):
    name = name.replace('\\', '/')
    bn = os.path.basename(name)
    x = bn.split('-')
    sDate = "%s-%s-%s" % (x[0], x[1], x[2])
    lDate.append(sDate)
    spc = wrap(readSpectrum(name))
    sp = spc.getProperties()
    pc = sp.getNumericProperty(epq.SpectrumProperties.FaradayBegin)
    lt = sp.getNumericProperty(epq.SpectrumProperties.LiveTime)
    sp.setDetector(det)
    spc.rename(sDate)
    if bDisplayAll:
        display(spc)
    comp = epq.Composition(epq.Element.Cu)
    sf = epq.SpectrumFitter8(det, comp, spc)
    rois = sf.getROIS()
    props = sp.getDetector().getCalibration().getProperties()
    offset = props.getNumericWithDefault(epq.SpectrumProperties.EnergyOffset, 0.0)
    gain   = props.getNumericWithDefault(epq.SpectrumProperties.EnergyScale, 10.0)
    coeffs = [ offset, gain]
    sf.setEnergyScale(epq.SpectrumFitter8.EnergyScaleFunction(coeffs, 2))
    sf.setResolution(epq.SpectrumFitter8.FanoNoiseWidth(6.0));
    sf.setMultiLineset(sf.buildWeighted(rois));
    results = sf.compute()
    # sf.setEnergyScale(epq.SpectrumFitter8.AltEnergyScaleFunction(coeffs))
    # print(dir(results))
    # ['FWHMatMnKa', '__class__', '__copy__', '__deepcopy__', '__delattr__',
    #  '__doc__', '__ensure_finalizer__', '__eq__', '__format__', 
    #  '__getattribute__', '__hash__', '__init__', '__ne__', '__new__', 
    #  '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__str__',
    #  '__subclasshook__', '__unicode__', 'bremsstrahlungSpectrum', 
    #  'bremstrahlungModel', 'calibration', 'channelWidth', 'class',
    #  'equals', 'fanoFactor', 'fitSpectrum', 'format3Col', 
    #  'getBremsstrahlungSpectrum', 'getCalibration', 'getChannelWidth',
    #  'getClass', 'getFWHMatMnKa', 'getFanoFactor', 'getFitSpectrum', 
    #  'getFitTransitions', 'getGaussianWidth', 'getIntegratedIntensity',
    #  'getNoise', 'getPosition', 'getResidual', 'getSubSpectrum', 
    #  'getTransitions', 'getUnknown', 'getZeroOffset', 'hashCode',
    #  'noise', 'notify', 'notifyAll', 'residual', 'resolutionToHTML',
    #  'setBremsstrahlungModel', 'setBremsstrahlungSpectrum',
    #  'setBremstrahlungModel', 'setEnergyCalibration', 'setResolution',
    #  'tabulateResults', 'toHTML', 'toString', 'toTransitionString',
    #  'transitions', 'unknown', 'wait', 'zeroOffset']
    for i in range(mIterations-1):
        results = sf.recompute(10.0, 0.3)

    # fit = results.fitSpectrum
    best = wrap(sf.getBestFit())
    best.rename('Best Fit-%s' % sDate)
    if bDisplayAll:
        display(best)
    cu = wrap(results.fitSpectrum)
    cu.rename('Cu-%s' % sDate)
    display(cu)

    cw = results.getChannelWidth()
    lCwMu.append(cw.floatValue())
    lCwUN.append(cw.uncertainty())

    # print(results.getChannelWidth().toString())
    zo = results.getZeroOffset()
    lZoMu.append(zo.floatValue())
    lZoUN.append(zo.uncertainty())
    
    res = results.FWHMatMnKa
    lRsMu.append(res.floatValue())
    lRsUN.append(res.uncertainty())
    
    # print(results.FWHMatMnKa.toString())
    intCuLa = epq.SpectrumUtils.backgroundCorrectedIntegral(cu, 740.0, 1040.0)
    lPkMu.append(intCuLa[0]/(pc*lt))
    lPkUN.append(intCuLa[1]/(pc*lt))

# write the data to a .csv file
# open file and write header
fi = open(csvFil, 'w')
strLine  = "date, " 
strLine += "cw.ev.mu, cw.ev.unc, "
strLine += "zo.ev.mu, zo.ev.unc, "
strLine += "mn.res.mu, mn.res.unc, "
strLine += "cu.la.cts.per.na.sec.mu, cu.la.cts.per.na.sec.unc\n"
fi.write(strLine)

#loop through the arrays
for i in range(len(lDate)):
    strLine  = "%s, " % (lDate[i])
    strLine += "%.6f, %.6f, " % (lCwMu[i], lCwUN[i])
    strLine += "%.5f, %.5f, " % (lZoMu[i], lZoUN[i])
    strLine += "%.4f, %.4f, " % (lRsMu[i], lRsUN[i])
    strLine += "%.2f, %.2f\n" % (lPkMu[i], lPkUN[i])
    fi.write(strLine)

fi.close()




# clean up cruft
shutil.rmtree(rptDir)
print "Done!"

end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg
