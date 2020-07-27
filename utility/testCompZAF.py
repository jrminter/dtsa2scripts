# -*- coding: utf-8 -*-
"""
testCompZAF.py

"""

import os, shutil

start = time.time()

det = findDetector("Oxford p4 05eV 2K")
e0  = 3.     # keV
lt  = 100     # sec
pc  =   5.0   # nA
withPoisson = True

stdName = "SiO2"
sio2 = material("SiO2", density=2.65)

unkName ="Si"
si  = material("Si", density=2.648)


"""
From the certificate for SRM 470

The values in the dictionary are mass fractions of the component
oxides.

std.name = 'K-412'
mix = mixture({ "SiO2"  : 0.4535,
                "MgO"   : 0.1933,
                "CaO"   : 0.1525,
                "FeO"   : 0.0996,
                "Al2O3" : 0.0927},
               density=2.600,
               std.name=std.name)
"""

# should not need to change below here

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
outDir = homDir + "/Documents/git/dtsa2Scripts/html"
rptDir = wrkDir + '/testCompZAF Results/'

dose = pc * lt  # na-sec"

DataManager.clearSpectrumList()

"""

std = simulate(sio2, det, e0, dose, True)
fmtS = "%s-%g-kV"
sName = fmtS % (sio2.name, e0)
std.rename(sName)
std.setAsStandard(sio2)
std.display()
sp = std.getProperties()

"""

sp = epq.SpectrumProperties()


"""
XPP1991 - Chantler2005

Algorithm   XPP - Pouchou & Pichoir Simplified
MAC NIST-Chantler 2005
E0  15 keV

IUPAC   Seigbahn    Standard    Energy   ZAF      Z   A   F k-ratio
Si K-L3 Si Kα1  SiO2    1.7397  1.1761  1.0515  1.1185  1.0000  2.516151
Si K-M3 Si Kβ1  SiO2    1.8290  1.1602  1.0515  1.1033  1.0000  2.482053
"""

print("XPP1991 - Chantler2005")
resXPP = zaf(si, det, e0, epq.XPP1991(),
             epq.MassAbsorptionCoefficient.Chantler2005,
             xtra=sp, stds={"Si": sio2}, mode="EDS")


"""
PAP1991-Chantler2005

Algorithm   PAP - Pouchou & Pichoir's Full φ(ρz)
MAC NIST-Chantler 2005
E0  15 keV

IUPAC   Seigbahn    Standard    Energy   ZAF      Z   A   F k-ratio
Si K-L3 Si Kα1  SiO2    1.7397  1.1782  1.0515  1.1204  1.0000  2.520465
Si K-M3 Si Kβ1  SiO2    1.8290  1.1617  1.0515  1.1048  1.0000  2.485316
"""

print("\n\n")

print("PAP1991-Chantler2005")
resXPP = zaf(si, det, e0, epq.PAP1991(),
             epq.MassAbsorptionCoefficient.Chantler2005,
             xtra=sp, stds={"Si": sio2}, mode="EDS")





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
