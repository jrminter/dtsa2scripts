# -*- coding: utf-8 -*-
"""
test-material-function.py

Author: J. R. Minter
Started: 2018-11-04
Modified:

Note: tabs set to 4 spaces

Density and composition for K412 from here:
https://probesoftware.com/smf/index.php?topic=108.msg391#msg391
Note: John Donovan uses oxide fractions

 ELEMENT   K-RAW K-VALUE ELEMWT% OXIDWT% ATOMIC% FORMULA KILOVOL
   Si ka  .00000  .15769  21.199  45.352  16.571    .412   15.00
   Fe ka  .00000  .06536   7.742   9.960   3.044    .076   15.00
   Mg ka  .00000  .07454  11.657  19.331  10.530    .262   15.00
   Ca ka  .00000  .10037  10.899  15.250   5.970    .149   15.00
   Al ka  .00000  .03226   4.906   9.270   3.992    .099   15.00
   Mn ka  .00000  .00064    .077    .099    .031    .001   15.00
   O                      43.597    .800  59.822   1.489
   Na ka  .00000  .00020    .043    .058    .041    .001   15.00
   TOTAL:                100.120 100.120 100.000   2.489


Done!
This script required 0.000 min
Elapse: 0:00:00.0

"""

import sys
import os
import time
import shutil
import gov.nist.microanalysis.Utility as epu
import gov.nist.microanalysis.EPQLibrary as epq
import dtsa2.mcSimulate3 as mc3
import dtsa2.jmGen as jmg
import java.util as ju

gitDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/utility"
prjDir = gitDir + relPrj
rptDir = prjDir + '/test-material-function Results/'
# spcDir = gitDir + '/dtsa2Scripts/spc/test-material-function'
# jmg.ensureDir(spcDir)

def create_material_from_mf(l_comp, l_mf, density, name):
    """
    create_material_from_mf(l_comp, l_mf, density, name)

    Create a material from lists of compositions and mass fractions
    and the density and a name.

    Input
    -----
    l_comp  A list of compositions.

    l_mf    A list of mass fractions.

    density A number. The density in g/cm3

    name    A string. The name for the material. Use letters and numbers
            without spaces or + or -. Simple mnemonics are better.

    Return
    ------
    A material

    Example
    -------
    # Note: composition came from J. Donovan:  Std 160 NBS K-412 mineral glass
    l_comps = [epq.Element.O,  epq.Element.Si, epq.Element.Mg, epq.Element.Ca,
               epq.Element.Fe, epq.Element.Al, epq.Element.Mn, epq.Element.Na]
    l_mfs   = [0.43597 ,0.21199 ,0.11657 ,0.10899,
               0.07742, 0.04906, 0.00077, 0.00043]
    mat = create_material_from_mf(l_comps, l_mfs, 2.66, "K412")
    print(mat)
    print(mat.toHTMLTable())
    """
    comp = epq.Composition(l_comp, l_mf)
    mat = epq.Material(comp, epq.ToSI.gPerCC(density))
    mat.setName(name)
    return(mat)


DataManager.clearSpectrumList()
start = time.time()

l_comps = [epq.Element.O,  epq.Element.Si, epq.Element.Mg, epq.Element.Ca,
           epq.Element.Fe, epq.Element.Al, epq.Element.Mn, epq.Element.Na]
l_mfs   = [0.43597 ,0.21199 ,0.11657 ,0.10899,
           0.07742, 0.04906, 0.00077, 0.00043]

mat = jmg.create_material_from_mf(l_comps, l_mfs, 2.66, "K412")
print(mat)
print(mat.toHTMLTable())
# print(dir(mat))

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




