# -*- coding: utf-8 -*-
"""
DTSA-II Script - J. R. Minter - 2016-10-12

testGetMassFractions.py

  Date      Who  Comment
----------  ---  -----------------------------------------------
2016-10-12  JRM  Mass fractions the easy way...

Elapse: 0:00:00.0  ROCPW7ZC5C42

"""

import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3,FluorescenceXRayGeneration3, XRayTransport3", None)

import os, shutil
import gov.nist.microanalysis.EPQLibrary as epq
import dtsa2.jmGen as jmg

gitDir  = os.environ['GIT_HOME']
relPrj  = "/dtsa2Scripts/utility"
prjDir  = gitDir + relPrj
rptDir  = prjDir + '/testGetMassFractions Results/'



def getMassFractions(compound, elemList, iDigits):
    """"
    getMassFractions(compound, elemList, iDigits)

    A utility function to compute the mass fractions for a compound

    Parameters
    ----------

    compound - string
        The stoichiometry as a molecular formula
    elemList = a list of type epq.element.symbol
        the elements in the compound
    iDigits - integer
        decimal places to round

    Returns
    -------
    A dictionary of {Symbol : mass-fraction}

    Example
    -------
    elements = [epq.Element.Al, epq.Element.Zn, epq.Element.O]
    massFra = getMassFractions("Al2Zn98O100", 5.61, elements, 5)
    """
    mat = material(compound)
    mf = {}
    for el in elemList:
        wf = round(mat.weightFractionU(el, True).doubleValue(), iDigits)
        mf[el.toAbbrev().encode('ascii','ignore')] = wf
    
    return(mf)

def summarizeMaterial(mat, iDigits=5):
    """"
    summarizeMaterial(mat, iDigits)

    A utility function to summarize a material

    Parameters
    ----------

    mat - a DTSA material
        The material to list
    iDigits - integer
        decimal places to round

    Returns
    -------
    A tuple: ( name, dictionary{Symbol : mass-fraction}, density)

    Example
    -------
    import gov.nist.microanalysis.EPQLibrary as epq
    kapton = epq.Material(epq.Composition([ epq.Element.C,
                                            epq.Element.O,
                                            epq.Element.N,
                                            epq.Element.H],
                                           [ 0.69113,
                                             0.20924,
                                             0.07327,
                                             0.02636 ]
                                          ),
                                          epq.ToSI.gPerCC(1.420))
    kapton.setName("Kapton")

    out = summarizeMaterial(kapton, 5)
    print(out)
    """
    density = round(mat.density/1000.0,iDigits) # in g/cc
    name = mat.name
    mf = {}
    elemList = mat.getElementSet()
    for el in elemList:
        wf = round(mat.weightFractionU(el, True).doubleValue(), iDigits)
        mf[el.toAbbrev().encode('ascii','ignore')] = wf

    rv = (name, mf, density)
    
    return(rv)

# elements = [epq.Element.Al, epq.Element.Zn, epq.Element.O]
# elements = [epq.Element.Zn, epq.Element.C, epq.Element.O, epq.Element.H]
# elements = [epq.Element.Al, epq.Element.N, epq.Element.O, epq.Element.H]

# test call from jmGen.py v 0.0.7
# strCmpd = "Al2Zn98O100"
# strCmpd = "ZnC4O6H10"

# AF2400
# elements = [epq.Element.C, epq.Element.O, epq.Element.F]
# strName = "AF2400"
# strCmpd = "C461O174F748"
# 
# AF1600
# elements = [epq.Element.C, epq.Element.O, epq.Element.F]
# strName = "AF1600"
# strCmpd = "C395O130F660"
#
# PTFE
# elements = [epq.Element.C, epq.Element.F]
# strName = "PTFE"
# strCmpd = "C2F4"

# elements = [epq.Element.Cu, epq.Element.O]
# strName = "Cu(OH)2"
# strCmpd = "Cu(OH)2"


# massFra = jmg.getMassFractions(strCmpd, elements, 5)

kapton = epq.Material(epq.Composition([ epq.Element.C,
                                        epq.Element.O,
                                        epq.Element.N,
                                        epq.Element.H],
                                       [ 0.69113,
                                         0.20924,
                                         0.07327,
                                         0.02636 ]
                                      ),
                                    epq.ToSI.gPerCC(1.420))
kapton.setName("Kapton")


# from Wikepedia from Corning and Schott
pyrex = epq.Material(epq.Composition([ epq.Element.B,
                                       epq.Element.O,
                                       epq.Element.Na,
                                       epq.Element.Mg,
                                       epq.Element.Al,
                                       epq.Element.Si,
                                       epq.Element.Cl,
                                       epq.Element.Ca,
                                       epq.Element.Fe],
                                       [ 0.03920,
                                         0.53850,
                                         0.03120,
                                         0.00030,
                                         0.01170,
                                         0.37720,
                                         0.00100,
                                         0.00070,
                                         0.00030]
                                      ),
                                    epq.ToSI.gPerCC(2.30))

pyrex.setName("Pyrex")



testOxide = material("0.5*SiO2+0.5*Fe2O3", density=3.200)
testOxide.setName("TestOxide")

# From Table 2 of SP260-74
k411 = material("0.5430*SiO2+0.1467*MgO+0.1547*CaO+0.1442*FeO", density=2.600)
k411.setName("K411")
out = jmg.summarizeMaterial(k411, 8)
print(out)

# (u'K411', 
#  'O': {'af': 0.60287639, 'wf': 0.42855296}, 
# 'Si': {'af': 0.20575277, 'wf': 0.25674404},
# 'Mg': {'af': 0.08286755, 'wf': 0.0894855},
# 'Ca': {'af': 0.06280718, 'wf': 0.11183761},
# 'Fe': {'af': 0.04569611, 'wf': 0.11337989}},
#  2.6)


k412 = material("0.4535*SiO2+0.1933*MgO+0.1525*CaO+0.0927*Al2O3+0.0996*FeO", density=2.600)
k412.setName("K412")
out = jmg.summarizeMaterial(k412, 8)

zns = material("ZnS", density=4.09)
zns.setName("ZnS")
out = jmg.summarizeMaterial(zns, 8)

sio2 = material("SiO2", density=2.65)
sio2.setName("SiO2")
out = jmg.summarizeMaterial(sio2, 8)
print(out)

sic = material("SiC", density=3.21)
sic.setName("SiC")
out = jmg.summarizeMaterial(sic, 8)
print(out)

mat = epq.Material(epq.Composition([epq.Element.Sn,
                                     epq.Element.Zn,
                                     epq.Element.S],
                                    [0.3946,
                                     0.3053,
                                     0.3028]),
                      epq.ToSI.gPerCC(4.6))
mat.setName("Sn_Zn_S")
out = jmg.summarizeMaterial(mat, 8)
print(out)

s310 = epq.Material(epq.Composition([epq.Element.C,
                                     epq.Element.Si,
                                     epq.Element.P,
                                     epq.Element.S,
                                     epq.Element.Cr,
                                     epq.Element.Mn,
                                     epq.Element.Fe,
                                     epq.Element.Ni],
                                    [0.0024,
                                     0.0150,
                                     0.0004,
                                     0.0003,
                                     0.2500,
                                     0.0200,
                                     0.5069,
                                     0.2050]),
                      epq.ToSI.gPerCC(7.89))
s310.setName("S310")
out = jmg.summarizeMaterial(s310, 8)
print(out)



# (u'K412', {
# El af           wf
# Al 0.04041414 0.04947715
#  O 0.59398097 0.43120208
# Si 0.16775488 0.21377747
# Mg 0.10659535 0.11755429
# Ca 0.06044228 0.10991362
# Fe 0.03081238 0.07807540



# print(strName, strCmpd, massFra)


# clean up cruft
shutil.rmtree(rptDir)
print "Done!"

