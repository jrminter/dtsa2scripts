"""
testPrz.py

"""

import os, shutil

start = time.time()

det = findDetector("Oxford p4 05eV 2K")
e0  = 15.     # keV
lt  = 100     # sec
pc  =   5.0   # nA
withPoisson = True

name = "ZnO"
zno = material("ZnO", density=5.606)


"""
From the certificate for SRM 470

The values in the dictionary are mass fractions of the component
oxides.

name = 'K-412'
mix = mixture({ "SiO2"  : 0.4535,
                "MgO"   : 0.1933,
                "CaO"   : 0.1525,
                "FeO"   : 0.0996,
                "Al2O3" : 0.0927},
               density=2.600,
               name=name)
"""

# should not need to change below here

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
outDir = homDir + "/Documents/git/dtsa2Scripts/html"
rptDir = wrkDir + '/testPrz Results/'

dose = pc * lt  # na-sec"

DataManager.clearSpectrumList()

out = phirhoz(zno, e0, det, rhoZmax=None, nSteps=100, alg=epq.PAP1991())


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
