# -*- coding: utf-8 -*-
"""
testMixtureOrthoclase.py

Define orthoclase mineral

"""

import os, shutil
import dtsa2.dt2Conv as dt2c

name = 'Orthoclase'
rho = 2.560
mix = mixture({ "K2O"    : 0.16925,
                "Al2O3"  : 0.18309,
                "SiO2"   : 0.64758
               },
               density=rho,
               name=name)


# should not need to change below here

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
outDir = homDir + "/Documents/git/dtsa2Scripts/html"
rptDir = wrkDir + '/testMixtureOrthoclase Results/'

fiPath = outDir + "/%s.html" % (name)


f=open(fiPath, 'w')
f.write(mix.toHTMLTable())
f.close()

dt2c.addMatToDatabase(mix, name, rho)

# clean up cruft
shutil.rmtree(rptDir)
print "Done!"

