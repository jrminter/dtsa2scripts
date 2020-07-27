# -*- coding: utf-8 -*-

# util-wrapper.py

import sys
import os
import time
import shutil
import gov.nist.microanalysis.Utility as nu

gitDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/utility"
prjDir = gitDir + relPrj
rptDir = prjDir + '/util-wrapper Results/'


start = time.time()



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



