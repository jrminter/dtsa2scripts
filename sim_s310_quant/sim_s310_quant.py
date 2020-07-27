# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

sim_s310_and_quant.py

A reproducible script to compute spectra and quantify S310 steel

  Date      Who   Ver    Comment
----------  ---  -----   ------------------------------------------
2020-06-10  JRM  0.0.1   Initial test
                         For 10,000 traj
                         This script required 56.903 min
                         Elapse: 0:56:54.1
                         Note the difference in the SRM-1155 spectrum

"""

__revision__ = "$Id: sim_310_steel_and_quant_spectra.py John R. Minter $"
__version__ = "0.0.1"

import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3,BremsstrahlungXRayGeneration3,FluorescenceXRayGeneration3, XRayTransport3", None)

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

DataManager.clearSpectrumList()

def defineMat(elms,qty,name,density=None):
	c=epq.Composition(map(element,elms),qty,name)
	if density:
		c=epq.Material(c,epq.ToSI.gPerCC(density))
	return(c)

e0      = 20.0
c_thick = 0.020 # 20 nm
n_traj  = 10000
dose    = 120.0
det     = findDetector("Oxford p4 05eV 4K")

"""
Certified Mass Fraction Values for SRM 1155 Stainless Steel (Cr 18-Ni 12-Mo 2)
Constituent Mass Fraction
(%)
Coverage Factor, k
Z  Mass Fraction
-- -----------------------
C   0.0445 +/- 0.0014 3.20
Co  0.1052 +/- 0.0057 2.36
Cr 18.37   +/- 0.21 4.30
Cu  0.1734 +/- 0.0075 2.26
Mn  1.619  +/- 0.075 3.20
Mo  2.386  +/- 0.024 2.23
Ni 12.35   +/- 0.22 4.30
S   0.0175 +/- 0.0032 4.30
V   0.0508 +/- 0.0034 2.36 


"""



homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/sim_s310_quant"

ctg     = material("C", 1.8)

ni2p    = material("Ni2P", 7.39)

fe      = material("Fe", 7.874)

s310    = defineMat(("P",   "S",   "Cr",  "Mn",  "Fe",  "Ni", ),
	                (0.0004,0.0003,0.2500,0.0200,0.5069,0.2050),
	                "S310", 7.89)

srm1155 = defineMat(("C"   , "Co",   "Cr",  "Cu",   "Mn" , "Mo",  "Ni",  "S",   "V",),
	                (0.0445, 0.1052, 18.37, 0.1734, 1.619, 2.386, 12.35, 0.0175, 0.0508),
	                "SRM-1155", 7.89)


# Fe
spc_fe = mc3.coatedSubstrate(ctg, c_thick, fe, det, e0,
                               True, n_traj, dose, True, True, {})
sName = "%g nm C on Fe at %g kV" % (1000*c_thick, e0)
spc_fe.rename(sName)
spc_fe.display()
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_fe.setAsStandard(fe)
spc_fe.save(fi)

# S310
spc_s310 = mc3.coatedSubstrate(ctg, c_thick, s310, det, e0,
                               True, n_traj, dose, True, True, {})
sName = "%g nm C on S310 at %g kV" % (1000*c_thick, e0)
spc_s310.rename(sName)
spc_s310.setAsStandard(s310)
spc_s310.display()
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_s310.save(fi)

# Ni2P
spc_ni2p = mc3.coatedSubstrate(ctg, c_thick, ni2p, det, e0,
                               True, n_traj, dose, True, True, {})
sName = "%g nm C on Ni2P at %g kV" % (1000*c_thick, e0)
spc_ni2p.rename(sName)
spc_ni2p.setAsStandard(ni2p)
spc_ni2p.display()
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_ni2p.save(fi)

#SRM-1155
spc_srm1155 = mc3.coatedSubstrate(ctg, c_thick, srm1155, det, e0,
                               True, n_traj, dose, True, True, {})
sName = "%g nm C on SRM-1155 at %g kV" % (1000*c_thick, e0)
spc_srm1155.rename(sName)
spc_srm1155.setAsStandard(ni2p)
spc_srm1155.display()
fi =  wrkDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (n_traj)
spc_srm1155.save(fi)

end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg





