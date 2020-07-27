# -*- coding: utf-8 -*-
# get_current_calibrations.pt
# John Minter 2020-07-24 John Minter

det1 = findDetector("Si(Li)")
print(det1)

det2 = findDetector("Oxford p4 05eV 2K")
print(det2)

det3 = findDetector("Oxford p4 05eV 4K")
print(det3)

det4 = findDetector("SDD (Medium, 4096)")
print(det4)

det5 = findDetector("Bruker 5eV")
print(det5)

det6 = findDetector("FEI CM20UT EDAX-RTEM")
print(det6)

det7 = findDetector("FEI FIB620 EDAX-RTEM")
print(det7)

"""
Results 2020-07-04

Probe
=====
Si(Li) - FWHM[Mn Kα]=132.0 eV - initial

Sirion
======
Oxford p4 05eV 2K - FWHM[Mn Kα]=130.2 eV - initial
Oxford p4 05eV 2K - FWHM[Mn Kα]=130.2 eV - initial
Oxford p4 05eV 4K - FWHM[Mn Kα]=130.1 eV - initial

NIST Mira3
==========
SDD (Medium, 4096) - FWHM[Mn Kα]=130.6 eV - initial

FEI CM20UT
==========
FEI CM20UT EDAX-RTEM - FWHM[Mn Kα]=140.7 eV - initial

FEI FIB 620
===========
FEI FIB620 EDAX-RTEM - FWHM[Mn Kα]=130.6 eV - initial

"""
