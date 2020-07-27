# sim-k496-phi-rho-z.py
#   Date      Who  Comment
# ----------  ---  -----------------------------------------------
# 2020-06-27  JRM  Simulate the phi-rho-z curve for K496 for
#                  e0 in 5, 10, 15, 20, 25, and 30kV
# 
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

gitHom = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts"
prjDir = gitHom + relPrj
wrkDir = prjDir + "/sim-k496"
jmg.ensureDir(wrkDir)
os.chdir(wrkDir)
przDir = wrkDir + '/prz/'
# Save spectra in sim directory
jmg.ensureDir(przDir)
# jmg.ensureDir(rptDir)


det     = findDetector("Oxford p4 05eV 4K")
e0      = 7      # kV
nSteps  = e0*100 # steps
rho     = 2.6      # sec

l_comps = [epq.Element.Al,  epq.Element.Mg, epq.Element.P, epq.Element.O]
l_mfs   = [ 0.06470,        0.06650,        0.32900,       0.5390]
k496    = jmg.create_material_from_mf(l_comps, l_mfs,  2.6, "K-496")
xrts    = mc3.suggestTransitions(k496)
a       = jmg.compPhiRhoZ(k496, det, e0, nSteps, xrts,    alg=epq.PAP1991(), base="pap-prz", outdir=przDir)


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