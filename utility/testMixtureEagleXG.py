# -*- coding: utf-8 -*-
"""
testMixtureEagleXG.py

Define Eagle-XG glass from midpoints mass fractions 
computed from mid-point mole fractions in patent.
"""

import os, shutil
import dtsa2.dt2Conv as dt2c

name = 'Eagle-XG'
rho = 2.380
mix = mixture({"SiO2"  : 0.6447825,
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
               density=rho,
               name=name)


# should not need to change below here

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
outDir = homDir + "/Documents/git/dtsa2Scripts/html"
rptDir = wrkDir + '/testMixtureEagleXG Results/'

fiPath = outDir + "/%s.html" % (name)


f=open(fiPath, 'w')
f.write(mix.toHTMLTable())
f.close()

dt2c.addMatToDatabase(mix, name, rho)

# clean up cruft
shutil.rmtree(rptDir)
print "Done!"

