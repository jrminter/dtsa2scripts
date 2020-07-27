# -*- coding: utf-8 -*-
"""
testMixtureAdmiraltyBrass.py

Define a mixture for Admiralty Brass and it to the database.
To permit replicating Goldstein2018a Fig. 17.37 p. 257
Simulated admiralty brass (69 % Cu, 30 % Zn, and 1 % Tin) for various
different selections of generation modes 

"""

import os, shutil
import dtsa2.dt2Conv as dt2c



"""
From:
text and here:
https://www.thoughtco.com/common-brass-alloys-and-their-uses-603706

densities
https://www.angstromsciences.com/density-elements-chart

The values in the dictionary are mass fractions of the component
oxides.
"""

name = 'Admiralty Brass'
rho = 0.69*8.86 + 0.30*7.13 + 0.01*7.31
mix = mixture({ "Cu"    : 0.69,
                "Zn"    : 0.30,
                "Sn"    : 0.01
               },
               density=rho,
               name=name)


# should not need to change below here

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
outDir = homDir + "/Documents/git/dtsa2Scripts/html"
rptDir = wrkDir + '/testMixtureAdmiraltyBrass Results/'

fiPath = outDir + "/%s.html" % (name)


f=open(fiPath, 'w')
f.write(mix.toHTMLTable())
f.close()

dt2c.addMatToDatabase(mix, name, rho)

# clean up cruft
shutil.rmtree(rptDir)
print "Done!"

