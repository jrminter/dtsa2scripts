# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

"""


__revision__ = "$Id: make_k412.py John R. Minter $"
__version__ = "0.0.1"

import sys
import os
import glob
import shutil
import fnmatch
import time
import math
import csv
import codecs

sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3, FluorescenceXRayGeneration3, XRayTransport3", None)
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.microanalysis.NISTMonte.Gen3 as nm3
import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQLibrary.Detector as epd
import gov.nist.microanalysis.Utility as epu
import gov.nist.microanalysis.EPQTools as ept
import java.io as jio
import java.lang as jl
import java.util as ju

from java.lang import Double


import dtsa2 as dt2
import gov.nist.microanalysis.dtsa2 as gdtsa2
import dtsa2.mcSimulate3 as mc3

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
    import dtsa2.jmGen as jmg
    # Note: composition came from J. Donovan:  Std 160 NBS K-412 mineral glass
    l_comps = [epq.Element.O,  epq.Element.Si, epq.Element.Mg, epq.Element.Ca,
               epq.Element.Fe, epq.Element.Al, epq.Element.Mn, epq.Element.Na]
    l_mfs   = [0.43597 ,0.21199 ,0.11657 ,0.10899,
               0.07742, 0.04906, 0.00077, 0.00043]
    mat = jmg.create_material_from_mf(l_comps, l_mfs, 2.66, "K412")
    print(mat)
    print(mat.toHTMLTable())
    """

import dtsa2.jmGen as jmg
# Note: composition came from J. Donovan:  Std 160 NBS K-412 mineral glass
l_comps = [epq.Element.O,  epq.Element.Si, epq.Element.Mg, epq.Element.Ca,
           epq.Element.Fe, epq.Element.Al, epq.Element.Mn, epq.Element.Na]
l_mfs   = [0.43597 ,0.21199 ,0.11657 ,0.10899,
           0.07742, 0.04906, 0.00077, 0.00043]
mat = jmg.create_material_from_mf(l_comps, l_mfs, 2.66, "K412")
print(mat)

print(mat.toHTMLTable())
