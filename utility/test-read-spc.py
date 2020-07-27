# -*- coding: utf-8 -*-

# test-read-spc.py

import sys
import os
import time
import shutil
import gov.nist.microanalysis.Utility as nu

gitDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/utility"
prjDir = gitDir + relPrj
rptDir = prjDir + '/test-read-spc Results/'

spcDir = prjDir + "/sim-quant-sio-w-sio2"


start = time.time()

fi = spcDir + "/SiO2 std.msa"
spc = readSpectrum(fi)
spc.display()

sp = spc.getProperties()
print(dir(sp))

comp = sp.getCompositionWithDefault(epq.SpectrumProperties.StandardComposition, None)
print(comp)



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



