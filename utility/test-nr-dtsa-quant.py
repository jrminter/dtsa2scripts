# -*- coding: utf-8 -*-

# Use NIST Demo Spec to test scripting DTSA-II quant
# This is adapted from N. Ritchie's presentation 
# titled "Scripting Composition using DTSA-II"
# Dated May 23, 2014
#
# Script date 2018-10-19  J. R. Minter
import sys
import os
import time
import shutil
import gov.nist.microanalysis.Utility as nu
import dtsa2.mcSimulate3 as mc3
import dtsa2.jmGen as jmg
# import gov.nist.microanalysis.EPQLibrary as epq

# directories
gitDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/utility"
prjDir = gitDir + relPrj
resDir = prjDir + '/nr-dtsa-quant-res'
jmg.ensureDir(resDir)
rptDir = prjDir + '/test-nr-dtsa-quant Results/'
spcDir = gitDir + '/dtsa2Scripts/spc/K411-K412'


DataManager.clearSpectrumList()
start = time.time()

det = findDetector("Bruker 5 eV")

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


mgo = material("MgO", density=3.6)
mgospc = load_and_prep_spectrum(spcDir, "MgO std.msa", mgo, det)


caf2 = material("CaF2", density=3.18)
caf2spc = load_and_prep_spectrum(spcDir, "CaF2 std.msa", caf2, det)

fe = material("Fe", density=7.87)
fespc = load_and_prep_spectrum(spcDir, "Fe std.msa", fe, det)


al2o3 = material("Al2O3", density=3.987)
al2o3spc = load_and_prep_spectrum(spcDir, "Al2O3 std.msa", al2o3, det)


k411 = material("K411", density= 2.946)
k411spc = load_and_prep_spectrum(spcDir, "K411 std.msa", k411, det)


k412 = material("K412", density= 2.946)
k412spc = load_and_prep_spectrum(spcDir, "K412 std.msa", k412, det)


res = quantify(k411spc,{"Si":sispc,
                        "Mg":mgospc,
                        "Ca":caf2spc,
                        "Fe":fespc,
                        "Al":al2o3spc,
                        "O":al2o3spc})
tmp = res.getComposition()
print(tmp.descriptiveString(False))

res = quantify(k412spc,{"Si":sispc,
                        "Mg":mgospc,
                        "Ca":caf2spc,
                        "Fe":fespc,
                        "Al":al2o3spc,
                        "O":al2o3spc})
tmp = res.getComposition()
print(tmp.descriptiveString(False))

# Note: I can print the measured composition but not the analyzed
# composition.

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



