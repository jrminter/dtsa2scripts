# Functions for computing compositional statistics
# Author: Nicholas W. M. Ritchie 
# Created: 14-May-2008
# Updated by JRM: 2018-10-22
# Note: tabs set to 4 spaces

import sys
import os
import time
import shutil
import gov.nist.microanalysis.Utility as epu
import dtsa2.jmGen as jmg
import java.util as ju

gitDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/utility"
prjDir = gitDir + relPrj
rptDir = prjDir + '/test-comp-stats Results/'
spcDir = gitDir + '/dtsa2Scripts/spc/K411-K412'
det = findDetector("Bruker 5 eV")

DataManager.clearSpectrumList()
start = time.time()


def compStats(comps, norm=0):
    """compStats(comps, norm=0)

    Compute compositional statistics for the list of Composition objects.

    Parameters
    ----------
    comps A list of Composition objects.
        Example: [comp1, comp2, comp3]
    norm An integer. Default: 0.
        If norm=1, the Composition objects are normalized to 100%; otherwise not.

   Returns
   -------
   None.  Prints results to the command line
   """
    elms=ju.TreeSet()
    for comp in comps:
        elms.addAll(comp.getElementSet())
    dss=dict()
    for elm in elms:
        dss[elm] = epu.DescriptiveStatistics()
    for comp in comps:
        for elm in elms:
            dss[elm].add(comp.weightPercent(elm,norm))
    return dss

def quantStats(specs, prop=epq.SpectrumProperties.MicroanalyticalComposition):
    """quantStatus(specs,prop)

    Computes the average Composition from the specified list of ScriptableSpectrum
    objects where the specified Composition property is extracted.

    Parameters
    ----------
    specs A list of spectra to quantify
          Example: [spec1, spec2,...,specN]

    props the property to quantify
          Default: epq.SpectrumProperties.MicroanalyticalComposition

    Returns
    -------
    None  Prints to the command line
    """
    comps=[]
    str = None
    for spec in specs:
        c=spec.getProperties().getCompositionWithDefault(prop,None)
        if c!=None:
            if str==None:
                str=spec.toString()
            else:
                str="%s, %s" % (str, spec.toString())
            comps=comps+[c]
    print "Quantitative statistics from %d spectra (%s)\n  [%s]" % (len(comps), prop, str)            
    cs=compStats(comps)
    print "Z\tElement\tAverage\tStdDev"
    print "\t\t (%)\t (%)"
    for elm, ds in cs.iteritems():
        print "%d\t%s\t%2.3f\t%2.3f" % ( elm.getAtomicNumber(), elm.toAbbrev(), ds.average()*100.0, ds.standardDeviation()*100.0 )

# define a function to laod and prepare and display a spectrum for analysis
def load_and_prep_spectrum(path, name, mat, det):
    fi = path + "/" + name
    spc = readSpectrum(fi)
    sp = spc.getProperties()
    sp.setDetector(det)
    spc.setAsStandard(mat)
    spc.display()
    return(spc)

si = material("Si", density=2.33)
sispc = load_and_prep_spectrum(spcDir, "Si std.msa", si, det)
fi = spcDir + "/Si std.xml"
sispc.toXML(fi)

mgo = material("MgO", density=3.6)
mgospc = load_and_prep_spectrum(spcDir, "MgO std.msa", mgo, det)
fi = spcDir + "/MgO std.xml"
mgospc.toXML(fi)


caf2 = material("CaF2", density=3.18)
caf2spc = load_and_prep_spectrum(spcDir, "CaF2 std.msa", caf2, det)
fi = spcDir + "/CaF2 std.xml"
caf2spc.toXML(fi)

fe = material("Fe", density=7.87)
fespc = load_and_prep_spectrum(spcDir, "Fe std.msa", fe, det)
fi = spcDir + "/Fe std.xml"
fespc.toXML(fi)


al2o3 = material("Al2O3", density=3.987)
al2o3spc = load_and_prep_spectrum(spcDir, "Al2O3 std.msa", al2o3, det)
fi = spcDir + "/Al2O3 std.xml"
al2o3spc.toXML(fi)


k411 = material("K411", density= 2.946)
k411spc = load_and_prep_spectrum(spcDir, "K411 std.msa", k411, det)
fi = spcDir + "/K411 std.xml"
k411spc.toXML(fi)


k412 = material("K412", density= 2.946)
k412spc = load_and_prep_spectrum(spcDir, "K412 std.msa", k412, det)
fi = spcDir + "/K412 std.xml"
k412spc.toXML(fi)

k411spc_1 = k411spc.subSample(20.0)
k411spc_2 = k411spc.subSample(20.0)
k411spc_3 = k411spc.subSample(20.0)

res411_1 = quantify(k411spc_1,{"Si":sispc,
                               "Mg":mgospc,
                               "Ca":caf2spc,
                               "Fe":fespc,
                               "Al":al2o3spc,
                               "O":al2o3spc})
# print("dir of qus obj")
#@ print(dir(res411_1 ))
tmp1 = res411_1.getComposition()
# print("dir of composition obj")
# print(dir(tmp1))
print(tmp1.descriptiveString(False))
k411spc_1.setAsMicroanalyticalComposition(tmp1)
k411spc_1.display()
fi = spcDir + "/K411 smp1.xml"
k411spc_1.toXML(fi)

res411_2 = quantify(k411spc_2,{"Si":sispc,
                               "Mg":mgospc,
                               "Ca":caf2spc,
                               "Fe":fespc,
                               "Al":al2o3spc,
                               "O":al2o3spc})
tmp2 = res411_2.getComposition()
print(tmp2.descriptiveString(False))
k411spc_2.setAsMicroanalyticalComposition(tmp2)
k411spc_2.display()
fi = spcDir + "/K411 smp2.xml"
k411spc_2.toXML(fi)


res411_3 = quantify(k411spc_3,{"Si":sispc,
                               "Mg":mgospc,
                               "Ca":caf2spc,
                               "Fe":fespc,
                               "Al":al2o3spc,
                               "O":al2o3spc})
tmp3 = res411_3.getComposition()
print(tmp3.descriptiveString(False))
k411spc_3.setAsMicroanalyticalComposition(tmp3)
k411spc_3.display()
fi = spcDir + "/K411 smp3.xml"
k411spc_3.toXML(fi)


quantStats([k411spc_1, k411spc_2, k411spc_3], prop=epq.SpectrumProperties.MicroanalyticalComposition)


k412spc_1 = k412spc.subSample(20.0)
k412spc_2 = k412spc.subSample(20.0)
k412spc_3 = k412spc.subSample(20.0)


res412_1 = quantify(k412spc_1,{"Si":sispc,
                               "Mg":mgospc,
                               "Ca":caf2spc,
                               "Fe":fespc,
                               "Al":al2o3spc,
                               "O":al2o3spc})
tmp1 = res412_1.getComposition()
print(tmp1.descriptiveString(False))
k412spc_1.setAsMicroanalyticalComposition(tmp1)
k412spc_1.display()
fi = spcDir + "/K412 smp1.xml"
k412spc_1.toXML(fi)


res412_2 = quantify(k412spc_2,{"Si":sispc,
                               "Mg":mgospc,
                               "Ca":caf2spc,
                               "Fe":fespc,
                               "Al":al2o3spc,
                               "O":al2o3spc})
tmp2 = res412_2.getComposition()
print(tmp2.descriptiveString(False))
k412spc_2.setAsMicroanalyticalComposition(tmp2)
k412spc_2.display()
fi = spcDir + "/K412 smp2.xml"
k412spc_2.toXML(fi)

res412_3 = quantify(k412spc_3,{"Si":sispc,
                               "Mg":mgospc,
                               "Ca":caf2spc,
                               "Fe":fespc,
                               "Al":al2o3spc,
                               "O":al2o3spc})
tmp3 = res412_3.getComposition()
print(tmp3.descriptiveString(False))
k412spc_2.setAsMicroanalyticalComposition(tmp3)
k412spc_3.display()
fi = spcDir + "/K412 smp3.xml"
k412spc_3.toXML(fi)

quantStats([k412spc_1, k412spc_2, k412spc_3], prop=epq.SpectrumProperties.MicroanalyticalComposition)








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
