"""
testMixture.py

"""

import os, shutil

name = 'K-412'

"""
From the certificate for SRM 470

The values in the dictionary are mass fractions of the component
oxides.
"""
mix = mixture({ "SiO2"  : 0.4535,
                "MgO"   : 0.1933,
                "CaO"   : 0.1525,
                "FeO"   : 0.0996,
                "Al2O3" : 0.0927},
               density=2.600,
               name=name)


# should not need to change below here

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
outDir = homDir + "/Documents/git/dtsa2Scripts/html"
rptDir = wrkDir + '/testMixture Results/'

fiPath = outDir + "/%s.html" % (name)

print(type(mix))

f=open(fiPath, 'w')
f.write(mix.toHTMLTable())
f.close()

# clean up cruft
shutil.rmtree(rptDir)
print "Done!"

