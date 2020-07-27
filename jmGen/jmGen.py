# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

jmGen.py

A series of wrapper scripts to make DTSA-II automation easy

Place this file in DTSA_ROOT/lib/dtsa2/
call with:

import dtsa2.jmGen as jmg

  Date      Who   Ver    Comment
----------  ---  -----   ------------------------------------------
2016-05-14  JRM  0.0.2   Updated for Jupiter
2016-07-18  JRM  0.0.3   Added function anaCNO
2016-07-19  JRM  0.0.4   Changed limits in anaCNO
2016-09-15  JRM  0.0.5   Set probe current and live time in
                         StripBackground and added anaSi and anCSi
2016-09-16  JRM  0.0.6   Added anaCCu
2016-10-12  JRM  0.0.7   Added getMassFractions
2017-01-23  JRM  0.0.8   Added addCompositionsToDatabase
2017-01-24  JRM  0.0.81  Edited help strings
2017-03-30  JRM  0.0.85  Added computePhiRhoReportZ to do a 
                         Phi-Rho-Z analysis and report in Z (um)
2017-04-11  JRM  0.0.86  Modify two prz functions to choose between
                         XPP119 and PAP1991 algorithms
2017-05-25  JRM  0.0.87  Added multiFilmBSE
2017-05-26  JRM  0.0.88  Figured out import problem
2017-06-24  JRM  0.0.89  Added summarizeMaterial
2017-07-05  JRM  0.0.90  summarizeMaterial returns both mass fraction 
                         and atom fraction
2017-07-06  JRM  0.0.91  Added a addMatToDatabase. This will overwrite
                         an existing (problem) entry...
2017-08-30  JRM  0.0.92  Added convenience functions writeMF and writeAF
2017-08-31  JRM  0.0.93  merged writeMF and writeAF into writeCompo
                         and added calcMAN.
2018-09-24  JRM  0.0.94  Modified compKRs to return a list of krs where
                         each kr has [mean, uncertainty]
2018-09-26  JRM  0.0.95  Added a function to print a title string for 
                         K-ratios and uncertainties output by the
                         compKRs function.
2018-10-16  JRM  0.0.96  Added pretty_print_unc_val function
2018-10-20  JRM  0.0.97  Added make_list_unc_val_string
2018-11-05  JRM  0.0.97  Added create_material_from_mf

"""

__revision__ = "$Id: jmGen.py John R. Minter $"
__version__ = "0.0.98"

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
    comp = epq.Composition(l_comp, l_mf)
    mat = epq.Material(comp, epq.ToSI.gPerCC(density))
    mat.setName(name)
    return(mat)

def make_list_unc_val_string(id, luv):
    """
    make_list_unc_val_string(id, luv)

    Make a formatted string from an ID string and a list of uncertain values.

    Input
    -----
    id  A number or a string that will be output as a string.

    luv A list of DTSA-II UncertainValue2 items. These will be printed
        as comma-delimited pairs with 6 digits following the decimal.

    Return
    ------
    A string with comma-delimited values with the ID and mean and uncertainty
    for each item in the list. This is suitable for writing output to a .csv
    file.

    Example:
    --------
    import dtsa2.jmGen as jmg
    import gov.nist.microanalysis.Utility as epu
    nmZnO1   = 40.1
    uvOKa1   = epu.UncertainValue2(0.269157,0.000126)
    uvZnLa1  = epu.UncertainValue2(0.259251,9.4e-05)
    uvSiKa1  = epu.UncertainValue2(0.654561,8.4e-05)
    l_uvals = [uvOKa1, uvZnLa1, uvSiKa1]
    out = make_list_unc_val_string(nmZnO1, l_uvals)
    print(out)
    
    1> 40.1, 0.269157, 0.000126, 0.259251, 0.000094, 0.654561, 0.000084
    """
    lv = len(luv)
    i = 0
    rv = "%s, " % (id)
    for uv in luv:
        rc = round(uv.doubleValue(), 6)
        uc = round(uv.uncertainty(), 6)
        
        if i == lv-1:
            rv += "%g, %.6f" % (rc, uc)
        else:
            rv += "%g, %.6f, " % (rc, uc)
        i += 1
    return(rv)

def pretty_print_unc_val(uv, n_digits):
    """
    pretty_print_unc_val(uv, n_digits)

    Round and print an uncertain value

    Input
    -----
    uv An uncertain value.
        The value to be formatted for printing
    n_digits An integer.
        The number of digits to round both the mean and uncertainty

    Return
    ------
    A string ready for printing

    Example
    -------
    import dtsa2.jmGen as jmg
    import gov.nist.microanalysis.Utility as nu
    val = nu.UncertainValue2(1.2345678, 0.0005678)
    print(jmg.pretty_print_unc_val(val, 5))
    """
    rc = round(uv.doubleValue(), n_digits)
    uc = round(uv.uncertainty(), n_digits)
    rv = u'%g \u00B1' % (rc)
    rv += ' %g' % (uc)
    return(rv)

def write_kr_titl_string(trs):
    """
    write_kr_titl_string(trs)

    Write a title string for the list of lists of k-ratios and
    uncertainties output from the compKRs function.

    Input
    -----
    trs The transition set used to compute the K-ratios

    Return
    ------
    A string with the title
    """
    titl = "["
    for tr in trs:
        trName = tr.toString().split("-")[0]
        titl += "["
        titl += trName + ".mu," + trName + ".unc],"
    titl += "]"
    titl = titl.replace("],]", "]]")
    titl = titl.replace(" ", "")
    titl = titl.replace("All", "K") # for light element when all there is is K
    titl = titl.replace(",", ", ")  # add spaces in title at the end
    return(titl)

def calcMAN(cf, name):
    """
    calcMAN(cf, name)

    Create a material object from a chemical formula and output
    the mean atomic number.

    Input
    -----
    cf: string
        The chemical formula (can also be a sum of oxides...)
    name: string
        name for the material

    Return
    ------
    none: prints the MAN to the command line

    Example
    -------
    import dtsa2.jmGen as jmg

    jmg.calcMAN("0.1933*MgO+0.0927*Al2O3+0.4535*SiO2+0.1525*CaO+0.0996*FeO", "K412")

    """
    mat = dt2.parseChemicalFormula(cf, 0, name)
    man = mat.meanAtomicNumber()
    line = "%s mean atomic number: %.3f" % (name, man)
    print(line)

def writeCompo(cf, density, name, norm=False):
    """
    writeCompo(cf, density, name, norm=False)

    Create a material from a string and print out the composition as
    mass fraction and atom fraction.

    Input
    -----
    cf: string
        The chemical formula (can also be a sum of oxides...)
    density: float
        The density [g/cm3]
    name: string
        name for the material
    norm: Boolean (False)
        Force normalization of mass fraction

    Return
    ------
    mat: DTSA material
        The material object

    Example
    -------
    import dtsa2.jmGen as jmg

    mat = jmg.writeCompo("0.1933*MgO+0.0927*Al2O3+0.4535*SiO2+0.1525*CaO+0.0996*FeO", 2.600, "K412", True)

    """

    mat = dt2.parseChemicalFormula(cf, density, name)
    if norm:
        mat.forceNormalization()
    se = mat.getSortedElements() # an array highest mf to lowest
    print("")
    print(name)
    print("element, mass-fraction, at-fraction")
    sumMF = 0.
    sumAF = 0.
    for el in se:
        line = "%s, " % (el)
        d0 = mat.weightFractionU(el, True)
        line += "%.6f, " % (round(d0.floatValue(), 6))
        sumMF += round(d0.floatValue(), 6)
        d1 = mat.atomicPercentU(el, True)
        line += "%.6f" % (round(d1.floatValue(), 6))
        sumAF += round(d1.floatValue(), 6)
        print(line)
    print("")
    line = "sum mass fractions = %.6f" % round(sumMF, 6)
    print(line)
    line = "sum atom fractions = %.6f" % round(sumAF, 6)
    print(line)
    print("")
    return mat

def addMatToDatabase(mat, name, density):
    """addMatToDatabase(mat, name, density)

    Add (or correct) an entry in the database.

    Parameters
    ----------

    mat - a DTSA material
        The material to add
    name - string
        The name to asssociate
    density - float
        the denisty in g/cm3

    Returns
    -------
    none. It just prints a message.

    Example
    -------
    import gov.nist.microanalysis.EPQLibrary as epq
    import dtsa2.jmGen as jmg
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

    jmg.addMatToDatabase(kapton, "Kapton", 1.4200)

    """
    mat.setName(name)
    mat.setDensity(epq.ToSI.gPerCC(density))
    dt2.Database.addStandard(mat)
    print("Added %s" % mat)

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
    A tuple: ( name, dictionary{Symbol : {mf : mass-fraction, af: atom-fraction}}, density)

    Example
    -------
    import gov.nist.microanalysis.EPQLibrary as epq
    import dtsa2.jmGen as jmg
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

    out = jmg.summarizeMaterial(kapton, 5)
    print(out)
    """
    density = round(mat.density/1000.0,iDigits) # in g/cc
    name = mat.name
    mf = {}
    elemList = mat.getElementSet()
    for el in elemList:
        wf = round(mat.weightFractionU(el, True).doubleValue(), iDigits)
        af = round(mat.atomicPercentU(el,True).doubleValue(), iDigits)
        mf[el.toAbbrev().encode('ascii','ignore')] = {"wf": wf, "af": af}

    rv = (name, mf, density)
    
    return(rv)

def multiFilmBSE(layers, e0=20.0, nTraj=1000, bOutFull=False,
                 outDir='C:/Temp/',fullOutFile='BSE.csv'):
    """multiFilmBSE(layers, det, e0=20.0, withPoisson=True,
                    nTraj=1000, dose=60.0, bOutFull=False,
                    outDir='C:/Temp/',fullOutFile='BSE.csv')

    Monte Carlo simulate the backscattered fraction from a multilayer
    thin film.

    Parameters
    ----------
    layers - an iterable list of [material,thickness].
        Note thematerials must have associated densities.
    e0 - float (20)
        The accekerating voltage in kV
    nTraj - integer
        The number of trajectories to run
    bOutFull - boolean (False)
        A flag to write the full data to a csv file
    outDir - string ('C:/Temp/')
        Path for full output
    fullOutFile - string ('BSE.csv')
        name for full output file


    Returns
    -------
    tuple - (layers, e0, nTraj, backscatterFraction, 
             forwardscatterFraction)

    Example
    -------
    import dtsa2.jmGen as jmg
    al      = material("Al",  density=2.70)
    c       = material("C",   density=2.267)
    zno     = material("ZnO", density=5.61)
    si      = material("Si",  density=2.3296)

    layers = [ [c,   20*1.0e-9],
               [zno, 20*1.0e-9],
               [al,  20*1.0e-9],
               [si,  50.0e-6]
             ]

    a = jmg.multiFilmBSE(layers, 7.0, nTraj=10000)
    print(a)

    """
    # myRP = dt2.reportPath()
    # myRP.replace("\\", "/")
    # print(myRP)
    monte = nm.MonteCarloSS()
    monte.setBeamEnergy(epq.ToSI.keV(e0))
    p = {}
    p["Layers"] = layers
    sr = monte.getChamber()
    origin = (0.0, 0.0, 0.0)
    pos = origin
    for (mat, thickness,) in layers:
        if thickness <= 0.0:
            raise "The layer thickness must be larger than zero."
        monte.addSubRegion(sr, mat, nm.MultiPlaneShape.createFilm([0.0, 0.0, -1.0], pos, thickness))
        pos = epu.Math2.plus(pos, [0.0, 0.0, thickness + 1.0e-12])
    bs0 = nm.BackscatterStats(monte, 100)
    monte.addActionListener(bs0)
    ann = nm.AnnularDetector(1.0e-3, 10, (0.0, 0.0, -1.0e-3), (0.0, 0.0, 1.0)) 
    monte.addActionListener(ann)
    monte.runMultipleTrajectories(nTraj)
    if bOutFull:
        fi = outDir + fullOutFile
        print(u"Results -> %s" % fi)
        fos = jio.FileOutputStream(fi)
        try:
            osw = jio.OutputStreamWriter(fos)
            osw.append("Parameters:\n")
            osw.append("E0\t%g keV\n" % e0)
            for k, v in p.iteritems():
                osw.append("%s\t%s\n" % (k, v))
            ann.dump(osw)
            osw.flush()
            bs0.dump(fos)
        finally:
            fos.close()
    return (layers, e0, nTraj, bs0.backscatterFraction(), bs0.forwardscatterFraction())




def addCompositionsToDatabase(compoList):
    """
    addCompositionsToDatabase(compoList)

    A utility function to add compositions to the DTSA-II database
    This is based on an exemplar from N. Ritchie.

    Parameters
    ----------
    compoList - List of type epq.Composition
        The list of compositions to be added to the database

    Returns
    -------
    None

    Example
    -------
    # Add Au/Pd and Ta coating materials to the database
    import dtsa2.jmGen as jmg
    aupd = epq.Material(epq.Composition([epq.Element.Au, epq.Element.Pd],
                                        [0.60,0.40]),
                                        epq.ToSI.gPerCC(0.6*19.30+0.4*11.9))
    aupd.setName("Au-Pd")
    ta = material("Ta", density=16.4)

    compoList = [aupd, ta]
    jmg.addCompositionsToDatabase(compoList)

    """
    for comp in compoList:
        ses = gdtsa2.DTSA2.getSession()
        if not ses.findStandard(comp.getName()):
            ses.addStandard(comp)

def getMassFractions(compound, elemList, iDigits):
    """
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
    import dtsa2.jmGen as jmg
    elements = [epq.Element.Al, epq.Element.Zn, epq.Element.O]
    massFra = jmg.getMassFractions("Al2Zn98O98", elements, 5)
    """
    mat = dt2.material(compound)
    mf = {}
    for el in elemList:
        wf = round(mat.weightFractionU(el, True).doubleValue(), iDigits)
        mf[el.toAbbrev().encode('ascii','ignore')] = wf
    
    return(mf)


def anaSi(spc, det, digits=2, display=True):
    """anaSi

    Strip continuum and analyze a spectrum for Si

    Parameters
    ----------
    spc: A DTSA-II scriptable spectrum
        The spectrum to process
    det: The detector

    Returns
    -------
    A dictionary of tuples for C and Si where each tuple is the
    peak integral (in counts/nA-sec)  and the uncertainty.
    """
    sp = spc.getProperties()
    pc = sp.getNumericProperty(epq.SpectrumProperties.FaradayBegin)
    lt = sp.getNumericProperty(epq.SpectrumProperties.LiveTime)
    dose = pc*lt
    sp.setDetector(det)
    bks = StripBackground(spc, det)
    if display:
        nam = spc.getProperties().getTextProperty(epq.SpectrumProperties.SpectrumDisplayName) + '-bks'
        bks.rename(nam)
        bks.display()
    siInt = compPeakIntegral(bks, 1.746, 240, 1)

    out = { "Si": (round(siInt[0]/dose, digits), round(siInt[1]/dose, digits)) }

    return out

def anaCSi(spc, det, digits=2, display=True):
    """anaCSi

    Strip continuum and analyze a spectrum for C and Si

    Parameters
    ----------
    spc: A DTSA-II scriptable spectrum
        The spectrum to process
    det: The detector

    Returns
    -------
    A dictionary of tuples for C and Si where each tuple is the
    peak integral (in counts/nA-sec)  and the uncertainty.
    """
    sp = spc.getProperties()
    pc = sp.getNumericProperty(epq.SpectrumProperties.FaradayBegin)
    lt = sp.getNumericProperty(epq.SpectrumProperties.LiveTime)
    dose = pc*lt
    sp.setDetector(det)
    bks = StripBackground(spc, det)
    if display:
        nam = spc.getProperties().getTextProperty(epq.SpectrumProperties.SpectrumDisplayName) + '-bks'
        bks.rename(nam)
        bks.display()
    cInt = compPeakIntegral( bks, 0.282, 180, 1)
    siInt = compPeakIntegral(bks, 1.746, 240, 1)

    out = { "C": (round(cInt[0]/dose, digits), round(cInt[1]/dose, digits)),
            "Si": (round(siInt[0]/dose, digits), round(siInt[1]/dose, digits))
          }

    return out

def anaCCu(spc, det, digits=2, display=True):
    """anaCCu

    Strip continuum and analyze a spectrum for C and Cu

    Parameters
    ----------
    spc: A DTSA-II scriptable spectrum
        The spectrum to process
    det: The detector

    Returns
    -------
    A dictionary of tuples for C and Cu where each tuple is the
    peak integral (in counts/nA-sec)  and the uncertainty.
    """
    sp = spc.getProperties()
    pc = sp.getNumericProperty(epq.SpectrumProperties.FaradayBegin)
    lt = sp.getNumericProperty(epq.SpectrumProperties.LiveTime)
    dose = pc*lt
    sp.setDetector(det)
    bks = StripBackground(spc, det)
    if display:
        nam = spc.getProperties().getTextProperty(epq.SpectrumProperties.SpectrumDisplayName) + '-bks'
        bks.rename(nam)
        bks.display()
    cInt = compPeakIntegral( bks, 0.282, 180, 1)
    cuInt = compPeakIntegral(bks, 0.917, 425, 1)

    out = { "C": (round(cInt[0]/dose, digits), round(cInt[1]/dose, digits)),
            "Cu": (round(cuInt[0]/dose, digits), round(cuInt[1]/dose, digits))
          }

    return out

def anaCNO(spc, det, digits=2, display=False):
    """anaCNO

    Strip continuum and analyze a spectrum for C, N, O

    Parameters
    ----------
    spc: A DTSA-II scriptable spectrum
        The spectrum to process
    det: The detector

    Returns
    -------
    A dictionary of tuples for C, N, and O where each tuple is the
    peak integral (in counts/nA-sec)  and the uncertainty.
    """
    sp = spc.getProperties()
    pc = sp.getNumericProperty(epq.SpectrumProperties.FaradayBegin)
    lt = sp.getNumericProperty(epq.SpectrumProperties.LiveTime)
    dose = pc*lt
    sp.setDetector(det)
    bks = StripBackground(spc, det)
    if display:
        nam = spc.getProperties().getTextProperty(epq.SpectrumProperties.SpectrumDisplayName) + '-bks'
        bks.rename(nam)
        bks.display()
    cInt = compPeakIntegral(bks, 0.282, 120, 1)
    nInt = compPeakIntegral(bks, 0.392, 100, 1)
    oInt = compPeakIntegral(bks, 0.525, 120, 1)

    out = { "C":(round(cInt[0]/dose, digits), round(cInt[1]/dose, digits)),
            "N":(round(nInt[0]/dose, digits), round(nInt[1]/dose, digits)),
            "O":(round(oInt[0]/dose, digits), round(oInt[1]/dose, digits))}

    return out


def fixSpcName(spc):
    """
    fixSpcName(spc)

    Remove the file extension from a spectrum name

    Parameters
    ----------
    spc: criptableSpectrum
        The spectrum to process

    Returns
    -------
    res: ScriptabelSpectrum
        The stripped spectrum
    """
    nam = spc.getProperties().getTextProperty(epq.SpectrumProperties.SpectrumDisplayName)
    print(nam)
    new = nam.split('.')[0]
    print(new)
    spc.rename(new)
    return spc

def StripBackground(spc, det=None):
    """
    StripBackground(spc, det=None)

    Strip the background from a spectrum

    Parameters
    ----------
    spc: ScriptableSpectrum
        The spectrum to process
    det: DTSA detector instance (None)
        If a detector is passed, the detector is set

    Returns
    -------
    res: ScriptabelSpectrum
        The stripped spectrum
    """
    sp = spc.getProperties()
    pc = sp.getNumericProperty(epq.SpectrumProperties.FaradayBegin)
    lt = sp.getNumericProperty(epq.SpectrumProperties.LiveTime)

    nam = spc.getProperties().getTextProperty(epq.SpectrumProperties.SpectrumDisplayName)
    new =  nam + '-bks'
    if det != None:
        spc.getProperties().setDetector(det)
    spc = epq.SpectrumUtils.applyZeroPeakDiscriminator(spc)
    spc = epq.PeakStripping.Clayton1987.getStrippedSpectrum(spc)
    spc = dt2.wrap(spc)
    sp = spc.getProperties()
    sp.setNumericProperty(epq.SpectrumProperties.FaradayBegin, pc)
    sp.setNumericProperty(epq.SpectrumProperties.FaradayEnd, pc)
    sp.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
    spc.rename(new)
    return spc

def computePhiRhoZ(mat, det, e0, xrts, csv, nSteps=200, usePAP=True):
    """
    computePhiRhoZ(mat, det, e0, xrts, csv, nSteps=200, usePAP=True):
    Compute a phi-rho-z curve

    Parameters
    ----------
    mat: a material
        A composition to compute
    det: a detector
        The detector to use
    e0: float
        The beam energy in kV
    xrts: XRayTransitionSet
        the transitions to compute
    csv: string
        Path to .csv file
    nSteps: int (200)
        The number of steps to compute
    usePAP: boolean (True)
        PAP1991 if true otherwise XPP1991

    Returns
    -------
    None - but writes output to file csv

    """
    rhoZmax = epq.ElectronRange.KanayaAndOkayama1972.compute(mat, epq.ToSI.keV(e0))
    if usePAP:
        alg = epq.PAP1991()
    else:
        alg = epq.XPP1991()
    sp = epq.SpectrumProperties(det.getProperties())
    sp.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
    res = "Idx,rhoz(mg/cm^2)"
    for xrt in xrts:
        res = "%s,G(%s),E(%s)" % (res, xrt, xrt)
    fi = open(csv, 'w')
    line = res + '\n'
    fi.write(line)

    for step in range(0, nSteps):
        rz = step * rhoZmax / nSteps
        res = "%d,%g" % (step, 100.0 * rz) # in mg/cm^2
        for xrt in xrts:
            alg.initialize(mat, xrt.getDestination(), sp)
            res = "%s,%g,%g" % (res, alg.computeCurve(rz), alg.computeAbsorbedCurve(xrt, rz))
        line = res+"\n"
        fi.write(line)
        # print(res)
    fi.close()

def computePhiRhoReportZ(mat, det, e0, xrts, csv, nSteps=200, usePAP=True):
    """
    computePhiRhoReportZ(mat, det, e0, xrts, csv, nSteps=200, usePAP=True)

    Compute a phi-rho-z curve and report a depth profile in Z
    The Z values are in microns.

    Parameters
    ----------
    mat: a material
        A composition to compute
    det: a detector
        The detector to use
    e0: float
        The beam energy in kV
    xrts: XRayTransitionSet
        the transitions to compute
    csv: string
        Path to .csv file
    nSteps: int (200)
        The number of steps to compute
    usePAP: boolean (True)
        PAP1991 if true otherwise XPP1991

    Returns
    -------
    None - but writes output to file csv

    """
    rhoZmax = epq.ElectronRange.KanayaAndOkayama1972.compute(mat, epq.ToSI.keV(e0))
    zMax = rhoZmax / mat.getDensity()
    if usePAP:
        alg = epq.PAP1991()
    else:
        alg = epq.XPP1991()
    
    sp = epq.SpectrumProperties(det.getProperties())
    sp.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
    res = "Z(um)"
    for xrt in xrts:
        res = "%s,G(%s),E(%s)" % (res, xrt, xrt)
    fi = open(csv, 'w')
    line = res + '\n'
    fi.write(line)

    for step in range(0, nSteps):
        rz = step * rhoZmax / nSteps
        dz = step * zMax / nSteps
        res = "%g" % (1.0e6 * dz) # in um
        for xrt in xrts:
            alg.initialize(mat, xrt.getDestination(), sp)
            res = "%s,%g,%g" % (res, alg.computeCurve(rz), alg.computeAbsorbedCurve(xrt, rz))
        line = res+"\n"
        fi.write(line)
        # print(res)
    fi.close()

def getSpectrumFromDataCube(hsVecCube, det, spcFile, x, y, mapTime, pc=1.0, bDebug=False):
    """getSpectrumFromDataCube(hsVecCube, det, spcFile, x, y, mapTime, pc=1.0, bDebug=False)
    Get a calibrated spectrum from an open hyperspectral vector datacube
    (hsVecCube) which was typically opened from a .rpl file. Calibrate the
    energy using DTSA detector (det) and the sum spectrum file (spcFile).
    Live time is computed from the mapTime (sec)
    Extract the spectrum from the coordinates x (width) and y (height).
    The probe current (pc) defaults to 1.0. A bDebug flag (default False)
    presents diagnostic messages.

    Return a DTSA-II scriptable spectrum."""
    ss = dt2.readSpectrum(spcFile)
    props=ss.getProperties()
    props.setDetector(det)
    cw = ss.getChannelWidth()
    zo = ss.getZeroOffset()
    nc = ss.getChannelCount()
    if(bDebug):
        strDetProps = "Detector properties: cw %.5f, zo %.3f" % (cw, zo)
        print(strDetProps)
    
    ss.getProperties().setNumericProperty(epq.SpectrumProperties.FaradayBegin, pc)
    ss.getProperties().setNumericProperty(epq.SpectrumProperties.FaradayEnd, pc)

    if(bDebug):
        display(ss)
    newSP = epq.SpectrumProperties()
    newSP.addAll(props);
    newSP.setTextProperty(epq.SpectrumProperties.SourceFile, "paint-sum-dtsa.msa")
    lt = newSP.getNumericWithDefault(epq.SpectrumProperties.LiveTime, Double.NaN)
    ### y,x in oxford map
    w = hsVecCube.getWidth()
    h = hsVecCube.getHeight()
    d = hsVecCube.getDepth()
    if ((x < w) and (y < h)):
        hsVecCube.seek(y,x) # Ba
        dat = hsVecCube.readDouble(d)
        sp = epq.SpectrumUtils.toSpectrum(cw, -zo, nc, dat)
        sp = dt2.wrap(sp)
        sp.setEnergyScale(newSP.getNumericProperty(epq.SpectrumProperties.EnergyOffset), newSP.getNumericProperty(epq.SpectrumProperties.EnergyScale))
        # props.setTimestampProperty(epq.SpectrumProperties.AcquisitionTime, ts)
        sp.setProbeCurrent(pc)
    
        lt = mapTime / (w*h)
        sp.setLiveTime(lt)
        newName = "location (%d,%d)" % (x,y)
        sp.rename(newName)
        dt2.display(sp)
    else:
        print("check x and y")
        sp = None
    return sp

def getSirionSiCpsPerNa(e0):
    """getSirionSiCpsPerNa(e0)
    Output the cps per nA for the Sirion Oxford EDS detector for a given e0.
    These values were determined for PT 6, 5 eV/ch, 2K channels
    in 2014-09-12-Cu-Si.
    Example:
    import dtsa2.jmGen as jmg
    a = jmg.getSirionSiCpsPerNa(7.0)"""
    val = 0.0
    if(e0 == 5.0):
        val = 1832.4
    if(e0 == 7.0):
        val = 3884.2
    if(e0 == 10.0):
        val = 7573.1
    if(e0 == 15.0):
        val = 13921.0
    if(e0 == 20.0):
        val = 20059.1
    if(e0 == 30.0):
        val = 29362.6
    return val
    
def getSirionCuCpsPerNa(e0):
    """getSirionCuCpsPerNa(e0)
    Output the cps per nA for the Sirion Oxford EDS detector for a given e0.
    These values were determined for PT 6, 5 eV/ch, 2K channels in 2014-09-12-QC
    Example:
    import dtsa2.jmGen as jmg
    a = jmg.getSirionCuCpsPerNa(7.0)"""
    val = 0.0
    if(e0 == 5.0):
        val = 2525.7
    if(e0 == 7.0):
        val = 3890.7
    if(e0 == 10.0):
        val = 5391.1
    if(e0 == 15.0):
        val = 6786.4
    if(e0 == 20.0):
        val = 7145.9
    if(e0 == 30.0):
        val = 6702.2
    return val

def reportComps(comps,path,fString="%s,%6.5f"):
    """reportComps(comps,path,,fString="%s,%6.5f")
    Outputs a collection of compositions (as weight fraction) to a .csv file.
    Example:
    import dtsa2.jmGen as jmg
    jmg.reportComps(comps,"C:/Temp/fooWf.csv")
    """
    f = open(path, 'w')
    all = set()
    for comp in comps:
        all = all | set(comp.getElementSet())
    tmp = "Material"
    for elm in all:
        tmp = "%s,%s" % (tmp, elm.toAbbrev())
    f.write(tmp+"\n")
    for comp in comps:
        tmp = str(comp)
        for elm in all:
            tmp = fString % (tmp, comp.weightFraction(elm, False))
        f.write(tmp+"\n")
    # close the file
    f.close()
    
def reportAtmPct(comps,path):
    """reportAtmPct(comps,path)
    Outputs a collection of compositions (as atm pct) to a .csv file.
    Example:
    import dtsa2.jmGen as jmg
    jmg.reportAtmPct(comps,"C:/Temp/fooWf.csv")
    """
    f = open(path, 'w')
    all = set()
    for comp in comps:
        all = all | set(comp.getElementSet())
    tmp = "Material"
    for elm in all:
        tmp = "%s,%s" % (tmp, elm.toAbbrev())
    f.write(tmp+"\n")
    for comp in comps:
        tmp = str(comp)
        for elm in all:
            tmp = "%s,%6.5f" % (tmp, comp.atomicPercent(elm))
        f.write(tmp+"\n")
    # close the file
    f.close()

def isNaN(num):
    """isNaN(num)
    Check if a number is NaN, returning True of False"""
    return num != num

def checkNaN(x):
    """checkNaN(x)
    This checks if a value (e.g. K-ratio) is NaN and sets the value to
    zero if it is. This really helps when writing data frames to be
    read by R."""
    if isNaN(x):
        x = 0.0
    return x
    
def ensureDir(d):
    """ensureDir(d)
    Check if the directory, d, exists, and if not create it."""
    if not os.path.exists(d):
        os.makedirs(d)
        
def clearAllSpectra():
    """clearAllSpectra()
    Clear all spectra from the data manager."""
    DataManager = dt2.DataManager.getInstance()
    DataManager.clearSpectrumList()
    
def compPeakIntegral(spc, e0, wid, digits=1):
    """compPeakIntegral( spc, e0, wid, digits=1)
    Compute the background-corrected peak interval for the line at e0 keV with an
    integration width of wid eV and return the results the desired (typically 1)
    following the decimal."""
    ev = 1000*e0
    hw = 0.5*wid
    start = ev - hw
    end = ev + hw
    res = epq.SpectrumUtils.backgroundCorrectedIntegral(spc, start, end)
    r0 = round(res[0],digits)
    r1 = round(res[1],digits)
    # note: under the hood, this returns the estimate, res[0] 
    # and the uncertainty, res[1]. Typically it is one digit.
    # r0    = '%.0f' % res[0]
    # r1    = '%.0f' % res[1]
    ret = [r0, r1]
    return ret
    
def updateCommonSpecProps(spc, det, name="", liveTime=-1, probeCur=-1, e0=-1, wrkDist=-1):
    """updateCommonSpecProps(spc, det, name="", liveTime=-1, probeCur=-1, e0=-1, wrkDist=-1)
    Update common spectrum properties
    spc - the spectrum (must supply)
    det - the detector (must supply)
    name - display name (will update if not empty)
    liveTime - acquisition live time (sec) (will update if > 0)
    probeCur - probe current (in nA) (will update if > 0)
    e0 - incident beam energy (in kV) (will update if > 0)
    wrkDist - working distance (in mm) (will update if > 0)
    This is helpful if they were not set properly in a GUI..."""
    spc = matchDet(spc, det)
    props=spc.getProperties()
    props.setDetector(det)
    l = len(name)
    if(l > 0) :
        props.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName, name)
    if(liveTime > 0) :
        props.setNumericProperty(epq.SpectrumProperties.LiveTime, liveTime)
    if(probeCur > 0) :
        props.setNumericProperty(epq.SpectrumProperties.FaradayBegin, probeCur)
        props.setNumericProperty(epq.SpectrumProperties.FaradayEnd, probeCur)
    if(e0 > 0) :
        props.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
    if(wrkDist > 0) :
        props.setNumericProperty(epq.SpectrumProperties.WorkingDistance, wrkDist)
    return spc

def compKRs(unSpc, stds, trs, det, e0, digits=5):
    """compKRs(unSpc, stds, trs, det, e0, digits=5)
    Performs a MLLSQ filter-fit of unknown spectrum (unSpc),
    to the standard spectra in a list of dictionaries of standards.
    Each standard in the list is a dictionary like
    {"El":element("C"), "Spc":cSpc} containing the element name
    and the standard spectrum. The function then computes the
    k-ratios for a list of transition sets (trs). The
    function returns a list of [k-ratio, uncertainty],
    rounded to the desired
    number of digits. The function performs an NaN test, setting
    NaNs to zero."""
    # Now set up the calc
    qa = epq.CompositionFromKRatios()
    ff=epq.FilterFit(det,epq.ToSI.keV(e0))
    l = len(stds)
    for i in range(0, l):
        st=stds[i]
        el=st["El"]
        sp=st["Spc"]
        ff.addReference(el,sp)
    krs=ff.getKRatios(unSpc)
    kr=[]
    n = len(trs)
    for i in range(0, n):
        tr=trs[i]
        k=round(checkNaN(krs.getKRatio(tr)), digits)
        # if NAN, returns zero
        if k == 0:
            # just set ku to 0
            ku = 0
        else:
            ku = round(krs.getKRatioU(tr).uncertainty(), digits)
        krvals=[]
        krvals.append(k)
        krvals.append(ku)
        kr.append(krvals)
    return kr

def cropSpec(spc, start=0, end=2048, bClear=True):
    """cropSpec(spc, start=0, end=2048, bClear=True)
    crop the spectrum (spc) starting with a starting and ending channel.
    This transfers the channel width, zero offset, and probe current
    required for microanalysis. Note: This only seems to work
    right for non-zero starts with no default detector, but not so sure now...
    
    Example:
    import dtsa2.jmGen as jmg
    niSpc = ept.SpectrumFile.open(niFile)[0]
    niSpc.getProperties().setNumericProperty(epq.SpectrumProperties.FaradayBegin,1.0)
    niSpc.getProperties().setNumericProperty(epq.SpectrumProperties.FaradayEnd,1.0)
    ws = wrap(niSpc)
    niSpc = jmg.cropSpec(ws, start=50, end=2048, bClear=True)
    """
    dt2.display(spc)
    props = spc.getProperties()
    nm = props.getTextProperty(epq.SpectrumProperties.SpectrumDisplayName)
    ts = props.getTimestampProperty(epq.SpectrumProperties.AcquisitionTime)
    cw = spc.getChannelWidth()
    zo = spc.getZeroOffset()
    lt = spc.liveTime()
    pc = spc.probeCurrent()
    cr = epq.SpectrumUtils.slice(spc, start, end-start)
    dt2.DataManager.removeSpectrum(spc)
    if bClear:
        dt2.clear()
    # note fix to get zero offset right
    sp = epq.SpectrumUtils.toSpectrum(cw, zo+start*cw, end-start, cr)
    rm = dt2.wrap(sp)
    # nzo = zo+start
    # rm = dt2.wrap(epq.SpectrumUtils.remap(sp, nzo, cw))
    # remap never worked right...
    props = rm.getProperties()
    props.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName, nm)
    props.setTimestampProperty(epq.SpectrumProperties.AcquisitionTime, ts)
    rm.setProbeCurrent(pc)
    rm.setLiveTime(lt)
    # dt2.display(sp)
    return rm

def clipSpec(spc, lcChan=0, hcChan=2048, bClear=True):
    """clipSpec(spc, lcChan=0, hcChan=2048, bClear=True)
    clip the spectrum (spc) starting with a starting and ending channel.
    This transfers the channel width, zero offset, and probe current
    required for microanalysis.
    
    Example:
    import dtsa2.jmGen as jmg
    niSpc = ept.SpectrumFile.open(niFile)[0]
    niSpc.getProperties().setNumericProperty(epq.SpectrumProperties.FaradayBegin,1.0)
    niSpc.getProperties().setNumericProperty(epq.SpectrumProperties.FaradayEnd,1.0)
    ws = wrap(niSpc)
    niSpc = jmg.clipSpec(ws, lcChan=20, hcChan=2048)
    """
    dt2.display(spc)
    props = spc.getProperties()
    nc = spc.getChannelCount()
    nm = props.getTextProperty(epq.SpectrumProperties.SpectrumDisplayName)
    ts = props.getTimestampProperty(epq.SpectrumProperties.AcquisitionTime)
    cw = spc.getChannelWidth()
    zo = spc.getZeroOffset()
    lt = spc.liveTime()
    pc = spc.probeCurrent()
    cr = epq.SpectrumUtils.slice(spc, 0, hcChan)
    # zero uneeded area
    for i in range(lcChan):
        cr[i] = 0.
    # for i in range(hcChan, nc):
    #    cr[i] = 0.
    if bClear:
        dt2.clear()
    dt2.DataManager.removeSpectrum(spc)
    # note fix to get zero offset right
    sp = epq.SpectrumUtils.toSpectrum(cw, 0, hcChan, cr)
    sp = dt2.wrap(sp)
    props = sp.getProperties()
    props.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName, nm)
    props.setTimestampProperty(epq.SpectrumProperties.AcquisitionTime, ts)
    sp.setProbeCurrent(pc)
    sp.setLiveTime(lt)
    # dt2.display(sp)
    return sp

def getKalphaEnergy(elmName):
    """getKalphaEnergy("Cu")
    Returns the transition energy (in keV) for the element's K-alpha
    line as the weighted sum of Ka1 and Ka2."""
    elm = dt2.element(elmName)
    xrt = epq.XRayTransition.KA1
    try:
        if epq.XRayTransition.exists(elm, xrt):
            tr = epq.XRayTransition(elm, xrt)
            en1 = epq.FromSI.keV(tr.getEnergy())
            wt1 = tr.getWeight(epq.XRayTransition.NormalizeKLM)
            xrt = epq.XRayTransition.KA2
            tr = epq.XRayTransition(elm, xrt)
            en2 = epq.FromSI.keV(tr.getEnergy())
            wt2 = tr.getWeight(epq.XRayTransition.NormalizeKLM)
            en = (wt1*en1+wt2*en2)/(wt1+wt2)
            return round(en, 4) 
        else:
            print "invalid element"
            return -1.0
    except jl.IllegalArgumentException:
        pass

def getLalphaEnergy(elmName):
    """getLalphaEnergy("Cu")
    Returns the transition energy (in keV) for the element's L-alpha line as
    the weighted sum of La1 and La2. Loaded from jrmFunctions.py."""
    elm = dt2.element(elmName)
    xrt = epq.XRayTransition.LA1
    try:
        if epq.XRayTransition.exists(elm, xrt):
            tr = epq.XRayTransition(elm, xrt)
            en1 = epq.FromSI.keV(tr.getEnergy())
            wt1 = tr.getWeight(epq.XRayTransition.NormalizeKLM)
            xrt = epq.XRayTransition.LA2
            tr = epq.XRayTransition(elm, xrt)
            en2 = epq.FromSI.keV(tr.getEnergy())
            wt2 = tr.getWeight(epq.XRayTransition.NormalizeKLM)
            en = (wt1*en1+wt2*en2)/(wt1+wt2)
            return round(en, 4) 
        else:
            print "invalid element"
            return -1.0
    except jl.IllegalArgumentException:
        pass

def getKLenergy(elmName):
    """getKLenergy(elmName)
    returns the weighted energy (in keV) for the element's K-alpha and
    L-alpha lines as a dictionary. Loaded from jrmFunctions.py."""
    enK = getKalphaEnergy(elmName)
    enL = getLalphaEnergy(elmName)
    res = {"K-alpha" : enK, "L-alpha" : enL}
    return res

def simBulkSpcCor(name, matl, det, e0=15, nTraj=100, lt=100, pc=1, noise=True):
    """simBulkSpcCor(name, matl, det, e0=15, nTraj=100, lt=100, pc=1, noise=True)
    Use the mc3 Monte Carlo functions to simulate bulk spectrum from a material with a live time lt
    and a probe current pc at a voltage e0. This does corrects for continuum fluorescence.
    Example:
    import dtsa2.jmGen as jmg
    det=findDetector("FEI FIB620 EDAX-RTEM")
    ni = material("Ni", density=8.90)
    niSpc = jmg.simBulkSpcCor("Ni-sim", ni, det, e0=15,nTraj=100, lt=60, pc=1)
    display(niSpc)
    """
    
    dose = lt*pc
    spc = dt2.wrap(mc3.simulate(matl, det, e0, True, nTraj, dose, True, True))
    props = spc.getProperties()
    sName = "%s-%gkV" % (name, e0)

    props.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName, sName)
    props.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
    props.setNumericProperty(epq.SpectrumProperties.FaradayBegin,pc)
    props.setNumericProperty(epq.SpectrumProperties.FaradayEnd,pc)
    props.setNumericProperty(epq.SpectrumProperties.BeamEnergy,e0)
    return spc

def simBulkSpcUnCor(name, matl, det, e0=15, nTraj=100, lt=100, pc=1, noise=True):
    """simBulkSpcUnCor(name, matl, det, e0=15, nTraj=100, lt=100, pc=1, noise=True)
    Use the standard Monte Carlo functions to simulate bulk spectrum from a material with a live time lt
    and a probe current pc at a voltage e0.
    This does not correct for continuum fluorescence.
    Example:
    import dtsa2.jmGen as jmg
    det=findDetector("FEI FIB620 EDAX-RTEM")
    ni = material("Ni", density=8.90)
    niSpc = jmg.simBulkSpcUnCor("Ni-sim", ni, det, e0=15,nTraj=100, lt=60, pc=1)
    display(niSpc)
    """
    origin=epu.Math2.multiply(1.0e-3, epq.SpectrumUtils.getSamplePosition(det.getProperties()))
    # create a simulator and initialize initialize it.

    monte=nm.MonteCarloSS()
    monte.setBeamEnergy(epq.ToSI.keV(e0))

    # create the substrate
    monte.addSubRegion(monte.getChamber(),matl,nm.MultiPlaneShape.createSubstrate([ 0.0, 0.0, -1.0],origin))

    # add event listener to model characteristic radiation
    xrel=nm.XRayEventListener2(monte,det)
    monte.addActionListener(xrel)
    brem=nm.BremsstrahlungEventListener(monte,det)
    monte.addActionListener(brem)
    
    det.reset()
    monte.runMultipleTrajectories(nTraj)
    # Get the spectrum and properties
    dose = lt*pc
    spc=det.getSpectrum(dose*1.0e-9 / (nTraj * epq.PhysicalConstants.ElectronCharge))
    if noise:
        spc=epq.SpectrumUtils.addNoiseToSpectrum(spc,1.0)
    props=spc.getProperties()
    sName="%s-%gkV" % (name, e0)

    props.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName, sName)
    props.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
    props.setNumericProperty(epq.SpectrumProperties.FaradayBegin,pc)
    props.setNumericProperty(epq.SpectrumProperties.FaradayEnd,pc)
    props.setNumericProperty(epq.SpectrumProperties.BeamEnergy,e0)
    return spc

def spcTopHatFilter(spc, det, e0, fw=150, norm=False):
    """spcTopHatFilter(spc, det, e0, fw=150, norm=False)
    Compute a top hat filter for spectrum spc with the given detector
    the specified kV (e0) with the desired filter width and normalization.
    """
    spc = matchDet(spc, det)
    rawName=spc.getProperties().getTextProperty(epq.SpectrumProperties.SpectrumDisplayName)
    fltName=rawName+"-thf"
    spc.getProperties().setNumericProperty(epq.SpectrumProperties.FaradayBegin, 1.0)
    spc.getProperties().setNumericProperty(epq.SpectrumProperties.FaradayEnd, 1.0)
    spc.getProperties().setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
    sw=dt2.wrap(spc)
    spc=sw
    thf=epq.FittingFilter.TopHatFilter(fw, det.getChannelWidth())
    fs=epq.FilteredSpectrum(spc, norm)
    fs.setFilter(thf)
    fs.getProperties().setTextProperty(epq.SpectrumProperties.SpectrumDisplayName, fltName)
    fsw=dt2.wrap(fs)
    return fsw
    
def compPhiRhoZ(comp, det, e0, nSteps=100, xrts=[], alg=epq.PAP1991(), base="pap-prz", outdir="./"):
    """compPhiRhoZ(comp, det, e0, nSteps=100, xrts=[], alg=epq.PAP1991(), base="pap-prz", outdir="./")
    Computes the ionization    as a function of dept for the composition
    (comp) with the specified detector (det) with the specified number
    of steps (nSteps). Algorithm choices are, the epq.XPP1991()
    (simplified Pouchou), epq.PAP1991() (full Pouchou and Pichoir)
    or epq.Proza96Base() (Bastin et al) algorithms. The results
    are written to a .csv file in the output directory (outdir)
    This reproduces the example in McSwiggen, 'Char. of sub-micrometer features with the FE-EPMA',
    EMAS2013, Fig 7 A and B p. 64
    Example:
    import dtsa2.jmGen as jmg
    e0  =  7
    nSteps = 500
    rho = 0.5*(7.874 + 8.908)
    feni = material("FeNi", density=rho)
    det = findDetector("FEI FIB620 EDAX-RTEM")
    trs=[epq.XRayTransition(epq.Element.Fe,    epq.XRayTransition.LA1), epq.XRayTransition(epq.Element.Ni,    epq.XRayTransition.LA1)]
    a = jmg.compPhiRhoZ(feni, det, e0, nSteps, xrts=trs,    alg=epq.PAP1991(), base="pap-prz", outdir="c:/temp/")
    """
    sName = comp.getName()
    rho = comp.getDensity() # kg/m3
    rho /= 1000.0 # g/cm3
    # print(rho)
    sFile = "%s-%s-%g-kv" % (sName, base, e0)
    print("Computing " + sFile)
    fName = outdir + sFile + ".csv"
    fi = open(fName,'w')
    sp = epq.SpectrumProperties(det.getProperties())
    sp.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
    if len(xrts) < 1:
        xrts = dt2.majorTransitions(comp, e0)
    
    rhoZmax = epq.ElectronRange.KanayaAndOkayama1972.compute(comp, epq.ToSI.keV(e0))
    # let's stick with DTSA-II default of writing results in mg/cm2.
    # Do the transform to z (nm) at plotting time
    res = "Idx,rhoz(mg/cm^2)"
    # res = "Idx, z(nm)"
    for xrt in xrts:
        res = "%s,G(%s),E(%s)" % (res, xrt, xrt)
    res = res + "\n"
    fi.write(res)
    for step in range(0, nSteps):
        rz = step * rhoZmax / nSteps
        res = "%d,%g" % (step, 100.0 * rz) # in mg/cm^2
        # res = "%d,%g" % (step, 1000000.0 * rz / rho) # in nm
        for xrt in xrts:
            alg.initialize(comp, xrt.getDestination(), sp)
            res = "%s,%g,%g" % (res, alg.computeCurve(rz), alg.computeAbsorbedCurve(xrt, rz))
        res=res+"\n"
        fi.write(res)
    fi.close()
    
def simMatlOnSubTEM(det, e0, matl, matlThick, subMat, subThick, lNames, nTraj=10000, lt=100, pc=1):
    """simMatlOnSubTEM(det, e0, matl, matlThick, subMat, subThick, lNames, nTraj=10000, lt=100, pc=1)
    Simulate a material on a thin substrate in the TEM for the detector (det),
    at the accelerating voltage (e0) keV, for the DTSA material (matl) with
    thickness (matlThick) in nm on the substrate DTSA material (subMat) with
    thickness (subThick) in nm with a list of names lNames = [matName, subName] by computing
    nTraj trajectories assuming a live time (lt, sec) and a probe current (pc, in nA) and
    simulating a noisy spectrum.
    Example:
    import dtsa2.jmGen as jmg
    det=findDetector("FEI CM20UT EDAX-RTEM")
    # create materials
    ago=epq.Material(epq.Composition([epq.Element.Ag, epq.Element.O],[0.930958,0.069042]), epq.ToSI.gPerCC(7.14))
    c=epq.MaterialFactory.createPureElement(epq.Element.C)
    spc = jmg.simMatlOnSubTEM(det, 200.0, ago, 100.0, c, 50.0, ["Ag2O","C"],    nTraj=10000, lt=100, pc=1)
    display(spc)
    """
    dose = lt*pc
    layTh=matlThick*1.0e-9
    subTh=subThick*1.0e-9
    # place sample at optimal location for the detector
    origin=epu.Math2.multiply(1.0e-3, epq.SpectrumUtils.getSamplePosition(det.getProperties()))
    # create a simulator and initialize it.
    monte=nm.MonteCarloSS()
    monte.setBeamEnergy(epq.ToSI.keV(e0))
    # create the film
    lay=nm.MultiPlaneShape.createFilm([0.0, 0.0, -1.0],[0.0, 0.0, 0.0], layTh)
    # create the substrate
    sub=nm.MultiPlaneShape.createFilm([0.0, 0.0, -1.0],[0.0, 0.0,layTh], subTh)
    monte.addSubRegion(monte.getChamber(), matl, lay)
    monte.addSubRegion(monte.getChamber(), subMat, sub)
    # add event listener to model characteristic radiation
    xrel=nm.XRayEventListener2(monte,det)
    monte.addActionListener(xrel)
    # ei=nm.EmissionImage.watchDefaultTransitions(xrel, 512, 2*layThick, origin)
    # add event listener to model bremsstrahlung
    brem=nm.BremsstrahlungEventListener(monte,det)
    monte.addActionListener(brem)
    # reset the detector and run the trajectories
    det.reset()
    monte.runMultipleTrajectories(nTraj)
    # Get the spectrum a    properties
    spec=det.getSpectrum(dose*1.0e-9 / (nTraj * epq.PhysicalConstants.ElectronCharge))
    noisy=dt2.wrap(epq.SpectrumUtils.addNoiseToSpectrum(spec,1.0))
    props=noisy.getProperties()
    spName = "Sim %g nm %s on %g nm %s" % (matlThick, lNames[0], subThick,    lNames[1])
    props.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName, spName)
    props.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
    props.setNumericProperty(epq.SpectrumProperties.FaradayBegin,pc)
    props.setNumericProperty(epq.SpectrumProperties.FaradayEnd,pc)
    props.setNumericProperty(epq.SpectrumProperties.BeamEnergy,e0)
    return(noisy)
    
def simBrehmTEM(det, e0, matl, matlThick, subMat, subThick, lNames, nTraj=10000, lt=100, pc=1):
    """simBrehmTEM(det, e0, matl, matlThick, subMat, subThick, lNames, nTraj=10000, lt=100, pc=1)
    Simulate the bremsstrahlung continuum spectrum for the detector (det),
    at the accelerating voltage (e0) keV, for the DTSA material (matl) with
    thickness (matlThick) in nm on the substrate DTSA material (subMat) with
    thickness (subThick) in nm    a list of names lNames = [matName, subName] 
    by computing nTraj trajectories assuming a live time (lt, sec) and a probe
    current (pc, in nA) and if desired, (addNoise) simulating a noisy spectrum.
    
    An example:
    import dtsa2.jmGen as jmg
    det=findDetector("FEI CM20UT EDAX-RTEM")
    # create materials
    ago=epq.Material(epq.Composition([epq.Element.Ag, epq.Element.O],[0.930958,0.069042]), epq.ToSI.gPerCC(7.14))
    c=epq.MaterialFactory.createPureElement(epq.Element.C)
    brehm = jmg.simBrehmTEM(det, 200.0, ago, 200.0, c, 50.0, ["Ag2O","C"], nTraj=10000, lt=100, pc=1)
    display(brehm)
    """
    dose = lt*pc
    layTh=matlThick*1.0e-9
    subTh=subThick*1.0e-9
    # create a simulator and initialize initialize it.
    monte=nm.MonteCarloSS()
    monte.setBeamEnergy(epq.ToSI.keV(e0))
    # create the film
    lay=nm.MultiPlaneShape.createFilm([0.0, 0.0, -1.0],[0.0, 0.0, 0.0], layTh)
    # create the substrate
    sub=nm.MultiPlaneShape.createFilm([0.0, 0.0, -1.0],[0.0, 0.0,layTh], subTh)
    monte.addSubRegion(monte.getChamber(),matl,lay)
    monte.addSubRegion(monte.getChamber(),subMat,sub)
    # add event listener to model bremsstrahlung
    brem=nm.BremsstrahlungEventListener(monte,det)
    monte.addActionListener(brem)
    # reset the detector and run the trajectories
    det.reset()
    monte.runMultipleTrajectories(nTraj)
    spec=det.getSpectrum(dose*1.0e-9 / (nTraj * epq.PhysicalConstants.ElectronCharge))
    noisy=epq.SpectrumUtils.addNoiseToSpectrum(spec,1.0)
    props=noisy.getProperties()
    spName = "Sim brehmstrahlung %g nm %s on %g nm %s" % (matlThick, lNames[0], subThick,    lNames[1])
    props.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName, spName)
    props.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
    props.setNumericProperty(epq.SpectrumProperties.FaradayBegin, pc)
    props.setNumericProperty(epq.SpectrumProperties.FaradayEnd, pc)
    props.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
    return dt2.wrap(noisy)
    
def measProbeCurrentFromCu(thisCu, stdCu, det, e0, digits=4):
    """measProbeCurrentFromCu(thisCu, stdCu, det, e0, digits=4)
    Measure the probe current for a session spectrum (thisCu) from
    the Cu standard (a probe current proxy) relative to the reference
    spectrum (stdCu) spectrum used to record standards. This assumes
    a detector (det) and e0 kV. The results are rounded to the desired digits.
    We use the Cu lines.
    
    Returns a dictionary:
    pcMu - the mean relative probe current (mean K-ratio for optimum TS)
    pcSE - the standard error for the relative probe current
    optTS - the optimal transition set for this kV
    ffMet - the filter fit metric, close to 0 is good, close to 1 is bad.
    
    An example:
    import dtsa2.jmGen as jmg
    res = jmg.measProbeCurrentFromCu(unCuSpc, rfCuSpc, det, e0)
    """
    thisCu = matchDet(thisCu, det)
    stdCu = matchDet(stdCu, det)
    qa = epq.CompositionFromKRatios()
    ff = epq.FilterFit(det,epq.ToSI.keV(e0))
    ff.addReference(dt2.element("Cu"),stdCu)
    krs = ff.getKRatios(thisCu)
    opt = krs.optimalDatum(dt2.element("Cu"))
    optTS = opt.toString()
    ku=krs.getKRatioU(opt)
    pcMu = round(ku.doubleValue(), digits)
    pcSE = round(ku.uncertainty(), digits)
    ffMet = round(ff.getFitMetric(thisCu), digits)
    # return a dictionary
    ret = {"pcMu": pcMu, "pcSE": pcSE, "optTS": optTS, "ffMet": ffMet }
    return ret

    
def avgSpectra(dir, names, det, e0, wrkDist, pc=1, resName="Avg", debug=False):
    """avgSpectra(dir, names, det, e0, wrkDist, pc=1, resName="Avg", debug=False)
    Compute the average spectrum from a list (names) of file names, assuming the individuals
    were recorded using the detector (det) at e0 kV with a working distance (wrkDist) and a
    probe current (pc) with default 1. Return return the average spectrum
    with a display (resName). The (debug) flag prints file names if positive.
    
    Example:
    import dtsa2.jmGen as jmg
    findDetector("FEI FIB620 EDAX-RTEM")
    theAvg = jmg.avgSpectra('C:\Temp\', ['Cu-12-1.spc','Cu-12-2.spc','Cu-12-3.spc'], det, 12, 17.1, 1.0, resName="Cu Std 12 kV")
    display(theAvg)
    """
    if(debug):
        print(names)
    nSpec = len(names)
    factor = 1.0 / float(nSpec)
    sPath = dir+names[0]
    sum = dt2.wrap(ept.SpectrumFile.open(sPath)[0])
    props = sum.getProperties()
    # use the live time and acquisition time of the first spectrum for the average
    liveTime = props.getNumericProperty(epq.SpectrumProperties.LiveTime)
    realTime = props.getNumericProperty(epq.SpectrumProperties.RealTime)
    ts = props.getTimestampProperty(epq.SpectrumProperties.AcquisitionTime)

    for i in range(1, nSpec):
        sPath = dir+names[i]
        tmp = dt2.wrap(ept.SpectrumFile.open(sPath)[0])
        sum += tmp
    avg=factor*sum
    
    props = avg.getProperties()
    props.setDetector(det)
    props.setNumericProperty(epq.SpectrumProperties.FaradayBegin, pc)
    props.setNumericProperty(epq.SpectrumProperties.FaradayEnd, pc)
    props.setNumericProperty(epq.SpectrumProperties.LiveTime, liveTime)
    props.setNumericProperty(epq.SpectrumProperties.RealTime, realTime)
    props.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
    props.setNumericProperty(epq.SpectrumProperties.WorkingDistance, wrkDist)
    props.setTextProperty(epq.SpectrumProperties.SpecimenDesc, "Averaged " + resName)
    props.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName, resName)
    props.setTimestampProperty(epq.SpectrumProperties.AcquisitionTime, ts)
    return avg

def getSpcAcqTime(spc):
    """getSpcAcqTime(spc)
    Return a string representation of the the local acquisition date and time for the input spectrum.
    Example:
    import dtsa2.jmGen as jmg
    strTim = jmg.getSpcAcqTime(cuSpc)
    """
    props = spc.getProperties()
    ts = props.getTimestampProperty(epq.SpectrumProperties.AcquisitionTime)
    acqTim = ts.toLocaleString()
    return acqTim

def getSpcAcqTimeDT(spc):
    """getSpcAcqTimeDT(spc)
    Return a date/time object with the local acquisition date and time for the input spectrum.
    Example:
    import dtsa2.jmGen as jmg
    dt = jmg.getSpcAcqTimeDT(cuSpc)
    """
    props = spc.getProperties()
    ts = props.getTimestampProperty(epq.SpectrumProperties.AcquisitionTime)
    return ts

 
def measRefProbeCur(baseDir, vkV, det, rptFil, debug=False):
    """measRefProbeCur(baseDir, vkV, det, rptFil, debug=False)
    Measure the relative prove current for Cu reference spectra in the
    baseDir/spc/refs/%gkV/ directory using the Cu standard spectra in
    baseDir/msa/std/%gkV/Cu.msa
    for a list (vkV) of e0 values for the detector det.
    Write the report (a .csv file) to rptFil.
    
    Example:
    import dtsa2.jmGen as jmg
    jmg.measRefProbeCur(baseDir, vkV, det, rptFil)
    """
    # prepare the output file
    f=open(rptFil, 'w')
    strLine = 'e0, refFile, acqTime, pcMu, pcSE\n'
    f.write(strLine)

    for e0 in vkV:
        stdFil = baseDir + "/msa/std/%gkV/Cu.msa" % e0
        std = dt2.wrap(ept.SpectrumFile.open(stdFil)[0])
        std = matchDet(std, det)
        refDir = baseDir + "/spc/refs/%gkV/" % e0
        refQue = baseDir + "/spc/refs/%gkV/*.spc" % e0
        a = glob.glob(refQue)
        l = len(a)
        if ( l > 0):
            for i in range(l):
                ref = dt2.wrap(ept.SpectrumFile.open(a[i])[0])
                refName = os.path.basename(a[i])
                refTim = getSpcAcqTime(ref)
                ref = matchDet(ref, det)
                res = measProbeCurrentFromCu(ref, std, det, e0)
                if(debug):
                    print([res['pcMu'], res['pcSE']])
                strLine = "%s" % e0 + ", "
                strLine = strLine + "%s" % refName.replace(".spc", "") + ", "
                strLine = strLine + "%s" % refTim.replace(",", "") + ", "
                strLine = strLine + "%.5f" % res['pcMu'] + ", "
                strLine = strLine + "%.5f" % res['pcSE'] + "\n"
                f.write(strLine)

    f.close()

def makeAvgRefSpectra(prjBaseDir, refName, nDupl, vkV, det, wrkDist, pc=1, debug=False):
    """ makeAvgRefSpectra(prjBaseDir, refName, nDupl, vkV, det, wrkDist, pc=1, debug=False)
    Generate standard spectra for a project with base directory (prjBaseDir) for
    the standard with (refName) from (nDupl) duplicate spectra recorded at
    each of a list (vkV) of accelerating voltages with detector (det), working
    distance (wrkDist), and probe current (pc, default 1.0).
    A debug flag(default, False) prints names.
    This assumes EDAX spc spectra
    are all store in prjBaseDir with names of the form refName-12-1.spc where
    12 is the kV and 1 is the duplicate number. This writes the standards in
    files like prjBaseDir/msa/std/12kV/refName.msa.
    
    Example:
    import dtsa2.jmGen as jmg
    findDetector("FEI FIB620 EDAX-RTEM")
    vkV = [10, 12, 15, 20, 25, 30]
    prjBaseDir = 'C:/Temp'
    makeAvgRefSpectra(prjBaseDir, "Cu", 3, vkV, det, 17.1, 1.0, False)
    """
    for e0 in vkV:
        vNames = []
        disName='%s Ref %gkV' % (refName, e0)
        spcDir = prjBaseDir + '/spc/refs/'
        msaDir = prjBaseDir + '/msa/ref/%gkV/' % e0
        for i in range(nDupl):
            vNames.append('%s-%g-%d.spc' % (refName, e0, i+1))
        if(debug):
            print(vNames)
        
        theAvg = avgSpectra(spcDir, vNames, det, e0, wrkDist, pc, resName=disName, debug=debug)
        updateCommonSpecProps(theAvg, det, name=disName, probeCur=pc, e0=e0, wrkDist=wrkDist)
        dt2.display(theAvg)
        
        outFil = msaDir + '%s.msa' % refName
        a = glob.glob(outFil)
        if (len(a) > 0):
            os.remove(a[0])
        fos=jio.FileOutputStream(outFil)
        ept.WriteSpectrumAsEMSA1_0.write(theAvg,fos,ept.WriteSpectrumAsEMSA1_0.Mode.COMPATIBLE)
        fos.close()
        if(debug):
            print(outFil)

def measRefProbeCurMsa(baseDir, vkV, det, rptFil, debug=False):
    """measRefProbeCurMsa(baseDir, vkV, det, rptFil, debug=False)
    Measure the relative prove current for Cu reference spectra in the
    baseDir/msa/ref/%gkV/ directory using the Cu standard spectra in
    baseDir/msa/std/%gkV/Cu.msa
    for a list (vkV) of e0 values for the detector det.
    Write the report (a .csv file) to rptFil.
    
    Example:
    import dtsa2.jmGen as jmg
    jmg.measRefProbeCurMsa(baseDir, vkV, det, rptFil)
    """
    # prepare the output file
    f=open(rptFil, 'w')
    strLine = 'e0, refFile, acqTime, pcMu, pcSE\n'
    f.write(strLine)

    for e0 in vkV:
        stdFil = baseDir + "/msa/std/%gkV/Cu.msa" % e0
        std = dt2.wrap(ept.SpectrumFile.open(stdFil)[0])
        std = matchDet(std, det)
        refDir = baseDir + "/msa/ref/%gkV/" % e0
        refQue = baseDir + "/msa/ref/%gkV/*.msa" % e0
        a = glob.glob(refQue)
        l = len(a)
        if ( l > 0):
            for i in range(l):
                ref = dt2.wrap(ept.SpectrumFile.open(a[i])[0])
                refName = os.path.basename(a[i])
                refTim = getSpcAcqTime(ref)
                ref = matchDet(ref, det)
                res = measProbeCurrentFromCu(ref, std, det, e0)
                if(debug):
                    print([res['pcMu'], res['pcSE']])
                strLine = "%s" % e0 + ", "
                strLine = strLine + "%s" % refName.replace(".spc", "") + ", "
                strLine = strLine + "%s" % refTim.replace(",", "") + ", "
                strLine = strLine + "%.5f" % res['pcMu'] + ", "
                strLine = strLine + "%.5f" % res['pcSE'] + "\n"
                f.write(strLine)

    f.close()

def simAnaStdSpc(mat, e0, det, cuRef, magNoise=1.0):
    """simAnaStdSpc(mat, e0, det, cuRef, magNoise=1.0)
    Use an analytical model to simulate an analytical standard from the specified material(mat)
    at the accelerating voltage (e0 in kV) for the specified detector (det) given a Cu reference
    spectrum (cuRef) and adding noise (magNoise) with a magnitude that defaults to 1.0
    The function simulates a Cu spectrum, computes the correct probe current, and then
    simulates the standard.
    Example:
    import dtsa2.jmGen as jmg
    mat = material("Ag", density=10.49)
    det = findDetector("FEI FIB620 EDAX-RTEM")
    cuRef = dt2.wrap(ept.SpectrumFile.open("PathToRef")[0])
    spc = jmg.simAnaStdSpc(mat, 25, det, cuRef, magNoise=1.0)
    """
    cu = dt2.material("Cu", density=8.96)
    props = cuRef.getProperties()
    cuRef = matchDet(cuRef, det)
    # use the live time and acquisition time of the first spectrum for the average
    liveTime = props.getNumericProperty(epq.SpectrumProperties.LiveTime)
    ts = props.getTimestampProperty(epq.SpectrumProperties.AcquisitionTime)
    pc = props.getNumericProperty(epq.SpectrumProperties.FaradayBegin)
    
    sp = epq.SpectrumProperties(det.getProperties())
    sp.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
    sp.setNumericProperty(epq.SpectrumProperties.FaradayBegin, 1.0)
    sp.setNumericProperty(epq.SpectrumProperties.FaradayEnd, 1.0)
    sp.setNumericProperty(epq.SpectrumProperties.LiveTime, liveTime)
        
    res = dt2.wrap(epq.SpectrumSimulator.Basic.generateSpectrum(cu, sp, det, True))
    cuSim = dt2.wrap(epq.SpectrumUtils.addNoiseToSpectrum(res, 1.0))
    cuSim.rename("%s-Sim-%0.1f kV" % (cu, e0))
    
    a = measProbeCurrentFromCu(cuSim, cuRef, det, e0)
    pcMu = a['pcMu']
    factor = 1.0/pcMu
    print(factor)
    
    sp = epq.SpectrumProperties(det.getProperties())
    sp.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
    sp.setNumericProperty(epq.SpectrumProperties.FaradayBegin, pc)
    sp.setNumericProperty(epq.SpectrumProperties.FaradayEnd, pc)
    sp.setNumericProperty(epq.SpectrumProperties.LiveTime, liveTime)
    
    res = dt2.wrap(epq.SpectrumSimulator.Basic.generateSpectrum(mat, sp, det, True))
    res *= factor
    stdSim = dt2.wrap(epq.SpectrumUtils.addNoiseToSpectrum(res, magNoise))
    props = stdSim.getProperties()
    props.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
    props.setNumericProperty(epq.SpectrumProperties.FaradayBegin, pc)
    props.setNumericProperty(epq.SpectrumProperties.FaradayEnd, pc)
    props.setNumericProperty(epq.SpectrumProperties.LiveTime, liveTime)
    # set a time stamp for the time of the ref
    props.setTimestampProperty(epq.SpectrumProperties.AcquisitionTime, ts)
    stdSim.rename("%s" % mat)
    return stdSim

def simMcStdSpc(mat, e0, det, cuRef, nTraj=1000, debug=False):
    """simMcStdSpc(mat, e0, det, cuRef, nTraj=1000, debug=False)
    Use a Monte Carlo simulation with Brehmsstralung correction to simulate an 
    analytical standard from the specified material (mat)
    at the accelerating voltage (e0 in kV) for the specified detector (det) given a Cu reference
    spectrum (cuRef) and (nTraj) trajectories that defaults to 1.0. The debug
    flag (default False) will display intermediate results when True.
    The function simulates a Cu spectrum, computes the correct probe current, and then
    simulates the standard.
    Example:
    import dtsa2.jmGen as jmg
    mat = material("Ag", density=10.49)
    det = findDetector("FEI FIB620 EDAX-RTEM")
    cuRef = dt2.wrap(ept.SpectrumFile.open("PathToRef")[0])
    spc = jmg.simMcStdSpc(mat, 25, det, cuRef, nTraj=100)
    """
    cu = dt2.material("Cu", density=8.96)
    props = cuRef.getProperties()
    cuRef = matchDet(cuRef, det)
    # use the live time and acquisition time of the first spectrum for the average
    liveTime = props.getNumericProperty(epq.SpectrumProperties.LiveTime)
    ts = props.getTimestampProperty(epq.SpectrumProperties.AcquisitionTime)
    pc = props.getNumericProperty(epq.SpectrumProperties.FaradayBegin)
    
    sp = epq.SpectrumProperties(det.getProperties())
    sp.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
    sp.setNumericProperty(epq.SpectrumProperties.FaradayBegin, 1.0)
    sp.setNumericProperty(epq.SpectrumProperties.FaradayEnd, 1.0)
    sp.setNumericProperty(epq.SpectrumProperties.LiveTime, liveTime)
        
    cuSim = dt2.wrap(simBulkSpcCor("Cu-sim", cu, det, e0, nTraj, lt=liveTime, pc=1))
    if(debug):
        dt2.display(cuSim)
    
    a = measProbeCurrentFromCu(cuSim, cuRef, det, e0)
    pcMu = a['pcMu']
    factor = 1.0/pcMu
    print(factor)
    
    sp = epq.SpectrumProperties(det.getProperties())
    sp.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
    sp.setNumericProperty(epq.SpectrumProperties.FaradayBegin, pc)
    sp.setNumericProperty(epq.SpectrumProperties.FaradayEnd, pc)
    sp.setNumericProperty(epq.SpectrumProperties.LiveTime, liveTime)
    
    name = "%s" % mat
    stdSim = dt2.wrap(simBulkSpcCor( name, mat, det, e0, nTraj, lt=liveTime, pc=1))
    stdSim *= factor
    stdSim.rename(name)
    props = stdSim.getProperties()
    props.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
    props.setNumericProperty(epq.SpectrumProperties.FaradayBegin, pc)
    props.setNumericProperty(epq.SpectrumProperties.FaradayEnd, pc)
    props.setNumericProperty(epq.SpectrumProperties.LiveTime, liveTime)
    # set a time stamp for the time of the ref
    props.setTimestampProperty(epq.SpectrumProperties.AcquisitionTime, ts)
    if(debug):
        dt2.display(stdSim)
    return stdSim

def compareBulkSpc(unk, ref, det, e0, elem=epq.Element.Cu, trs="Cu K-L3"):
    """compareBulkSpc(unk, ref, det, e0, elem=epq.Element.Cu, trs="Cu K-L3")
    Compare a bulk spectrum (unk) to a reference (ref) from the same element (elem)
    recorded on the same detector (det) and kV (e0). Use a preferred transition
    set (trs). Useful for checking probe current changes.
    Example:
    import dtsa2.jmGen as jmg
    det=findDetector("FEI FIB620 EDAX-RTEM")
    res = jmg.compareBulkSpc(unk, ref, det, e0, elem=epq.Element.Cu, trs="Cu K-L3")
    """
    comp = epq.Composition(elem)
    ref = matchDet(ref, det)
    unk = matchDet(unk, det)
    sp = ref.getProperties()
    sp.setCompositionProperty(epq.SpectrumProperties.StandardComposition, comp)
    qus = epq.QuantifySpectrumUsingStandards(det, epq.ToSI.keV(e0))
    qus.addStandard(elem, comp, ref)
    xrt = dt2.transition(trs)
    for roi in qus.getRegionOfInterestSet(elem):
        if roi.contains(xrt):
            qus.setPreferredROI(elem, roi)
    a = qus.compute(unk)
    res = a.getComposition()
    wtFrMu = round(res.weightFraction(elem, False), 5)
    wtFrU    = res.weightFractionU(elem, False)
    out = [wtFrMu, wtFrU]
    return out

def readEdaxSpc(path, det, pc, wrkDist):
    """readEdaxSpc(path, det, pc, wrkDist)
    Read an EDAX .spc file from (path), set the DTSA detector (det), the
    probe current (pc) and working distance (wrkDist). This writes to a
    temp .msa file and re-reads to get consistency
    Example:
    import dtsa2.jmGen as jmg
    det = findDetector("FEI FIB620 EDAX-RTEM")
    spc = jmg.readEdaxSpc(path, det, 1.0, 17.2)
    """
    spc = dt2.wrap(ept.SpectrumFile.open(path)[0])
    dt2.display(spc)
    props = spc.getProperties()
    lt = props.getNumericProperty(epq.SpectrumProperties.LiveTime)
    e0 = props.getNumericProperty(epq.SpectrumProperties.BeamEnergy)
    name = props.getTextProperty(epq.SpectrumProperties.SpectrumDisplayName)
    props.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName, name)
    props.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
    props.setNumericProperty(epq.SpectrumProperties.FaradayBegin,pc)
    props.setNumericProperty(epq.SpectrumProperties.FaradayEnd,pc)
    props.setNumericProperty(epq.SpectrumProperties.BeamEnergy,e0)
    spc = matchDet(spc, det)
    props = spc.getProperties()
    props.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName, name)
    props.setDetector(det)
    props.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
    props.setNumericProperty(epq.SpectrumProperties.FaradayBegin, pc)
    props.setNumericProperty(epq.SpectrumProperties.FaradayEnd, pc)
    props.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
    return spc
    
def sumCounts(spc, start, end):
    """sumCounts(spc, start, end)
    Compute a raw integral for a spectrum by summing the counts from
    start to end channels. Returns the sum.
    Example:
    import dtsa2.jmGen as jmg
    a = jmg.sumCounts(s1, 0, 1500)
    print(a)
    """
    spc = dt2.wrap(spc)
    sum = 0.
    for i in range(start, end):
        cts = spc.getCounts(i)
        sum += cts
        # print sum     
    return sum    

def getCurrentDetectorCalibration(det):
    """getCurrentDetectorCalibration(det)
    Get the current calibration from the specified detector. Return a dictionary with
    the probe, name, det name, window, offset, scale, quad, fudge factor, solid angle, resolution, and guid
    Example:
    import dtsa2.jmGen as jmg
    det=findDetector("FEI FIB620 EDAX-RTEM")
    cal = jmg.getCurrentDetectorCalibration(det)
    print(cal)
    """
    prop = det.getProperties()
    cal    = det.getCalibration()
    prob = det.getOwner()
    name = det.getName()
    win    = det.getWindow()
    ad     = cal.getActiveDate().toString()
    off    = round(det.getZeroOffset(), 5)
    sca    = round(cal.getChannelWidth(), 5)
    quad = det.getQuadratic()
    fud    = cal.getFudgeFactor()
    soa    = round(prop.getNumericProperty(epq.SpectrumProperties.SolidAngle), 5)
    resn = round(prop.getNumericProperty(epq.SpectrumProperties.Resolution), 2)
    guid = prop.getTextProperty(epq.SpectrumProperties.DetectorGUID)
    res = {"Probe":prob, "Name": name, "Win":win, "ActiveDate": ad, "Offset":off, "Scale":sca,    "Quad":quad, "FudgeFac":fud, "sR": soa, "resn": resn, "GUID":guid}
    return res
    
def estimateProbeCurrentFromCu(cuSpc, det, iDigits=5):
    """estimateProbeCurrentFromCu(cuSpc, det)
    Estimate the probe current from a Cu standard spectrum recorded
    from a Cu standard from the current detector. The function simulates
    a bulk Cu spectrum from a Phi-Rho-Z model for the detector at the same
    e0 and live time, assuming a 1 nA probe current. The estimated probe
    current is the ratio of the total counts from the two spectra.
    Example:
    import dtsa2.jmGen as jmg
    det = findDetector("FEI FIB620 EDAX-RTEM")
    pc = jmg.estimateProbeCurrentFromCu(spc, det)
    """
    cuSpc = matchDet(cuSpc, det)
    lt = cuSpc.liveTime()
    cts = cuSpc.totalCounts()
    props = cuSpc.getProperties()
    e0 = props.getNumericProperty(epq.SpectrumProperties.BeamEnergy)
    # print(e0)
    cu = dt2.material("Cu", density=8.96)
    sim = dt2.simulate(cu, det=det, keV=e0, dose=lt, withPoisson=False)
    cSim = sim.totalCounts()
    uExp = math.sqrt(cts)/cts
    uSim = math.sqrt(cSim)/cSim
    pc = cts/cSim
    uTot = math.sqrt(uExp*uExp+uSim*uSim)*pc
    ret = {"nA":round(pc, iDigits), "sd":round(uTot,iDigits)}
    return ret

def makeStdSpcSpectra(prjBaseDir, stdName, rho, nDupl, vkV, vPc, det, wrkDist, debug=False):
    """ makeStdSpcSpectra(prjBaseDir, stdName, rho, nDupl, vkV, vPc, det, wrkDist, debug=False)
    Generate standard spectra for a project with base directory (prjBaseDir) for
    the standard with (stdName) and density (rho) g/cm3 from (nDupl) duplicate spectra recorded at
    each of a list (vkV) of accelerating voltages with a corresponding list of
    probe currents (vPc) with detector (det), working distance (wrkDist).
    A debug flag(default, False) prints names.
    This assumes EDAX .spc spectra that are all stored in prjBaseDir/msa/stds/ with
    names of the form stdName-12-1.spc where 12 is the kV and 1 is the duplicate number.
    This writes the standards in .msa format files like prjBaseDir/msa/std/12kV/stdName.msa.
    
    Example:
    import dtsa2.jmGen as jmg
    findDetector("FEI FIB620 EDAX-RTEM")
    vkV = [12] # , 15, 20, 25, 30]
    vPc = [0.1163] # , 0.1261, 0.1278, 0.1518, 0.1335]
    prjBaseDir = 'C:/Temp'
    jmg.makeStdSpcSpectra(prjBaseDir, "Cu", 8.96, 3, vkV, vPc, det, wrkDist, debug=False)
    """
    ensureDir(prjBaseDir + '/msa/')
    ensureDir(prjBaseDir + '/msa/std/')
    mat = dt2.material(stdName, density=rho)
    l = len(vkV)
    for i in range(l):
        e0 = vkV[i]
        pc = vPc[i]
        vNames = []
        disName='%s Std %gkV' % (stdName, e0)
        spcDir = prjBaseDir + '/spc/stds/'
        msaDir = prjBaseDir + '/msa/std/%gkV/' % e0
        ensureDir(msaDir)
        for i in range(nDupl):
            vNames.append('%s-%g-%d.spc' % (stdName, e0, i+1))
        if(debug):
            print(vNames)
        
        theAvg = avgSpectra(spcDir, vNames, det, e0, wrkDist, pc, resName=disName, debug=debug)
        updateCommonSpecProps(theAvg, det, name=disName, probeCur=pc, e0=e0, wrkDist=wrkDist)
        theAvg = matchDet(theAvg, det)
        theAvg.setAsStandard(mat)
        dt2.display(theAvg)
        
        outFil = msaDir + '%s.msa' % stdName
        a = glob.glob(outFil)
        if (len(a) > 0):
            os.remove(a[0])
        fos=jio.FileOutputStream(outFil)
        ept.WriteSpectrumAsEMSA1_0.write(theAvg,fos,ept.WriteSpectrumAsEMSA1_0.Mode.COMPATIBLE)
        fos.close()
        if(debug):
            print(outFil)

def makeStdMsaSpectra(prjBaseDir, stdName, rho, nDupl, vkV, vPc, det, wrkDist, debug=False):
    """ makeStdMsaSpectra(prjBaseDir, stdName, rho, nDupl, vkV, vPc, det, wrkDist, debug=False)
    Generate standard spectra for a project with base directory (prjBaseDir) for
    the standard with (stdName) and density (rho) g/cm3 from (nDupl) duplicate spectra recorded at
    each of a list (vkV) of accelerating voltages with a corresponding list of
    probe currents (vPc) with detector (det), working distance (wrkDist).
    A debug flag(default, False) prints names.
    This assumes MSA .msa spectra that are all stored in prjBaseDir/msa/stds/ with
    names of the form stdName-12-1.msa where 12 is the kV and 1 is the duplicate number.
    This writes the standards in files like prjBaseDir/std/12kV/stdName.msa.
    
    Example:
    import dtsa2.jmGen as jmg
    findDetector("FEI FIB620 EDAX-RTEM")
    vkV = [12, 15, 20, 25, 30]
    vPc = [0.1163, 0.1261, 0.1278, 0.1518, 0.1335]
    prjBaseDir = 'C:/Temp'
    jmg.makeStdMsaSpectra(prjBaseDir, "Ni", 8.90, 3, vkV, vPc, det, wrkDist, debug=False)
    """
    mat = dt2.material(stdName, density=rho)
    l = len(vkV)
    for i in range(l):
        e0 = vkV[i]
        pc = vPc[i]
        vNames = []
        disName='%s Std %gkV' % (stdName, e0)
        spcDir = prjBaseDir + '/msa/stds/'
        msaDir = prjBaseDir + '/std/%gkV/' % e0
        ensureDir(msaDir)
        for i in range(nDupl):
            vNames.append('%s-%g-%d.msa' % (stdName, e0, i+1))
        if(debug):
            print(vNames)
        
        theAvg = avgSpectra(spcDir, vNames, det, e0, wrkDist, pc, resName=disName, debug=debug)
        updateCommonSpecProps(theAvg, det, name=disName, probeCur=pc, e0=e0, wrkDist=wrkDist)
        theAvg.setAsStandard(mat)
        dt2.display(theAvg)
        
        outFil = msaDir + '%s.msa' % stdName
        a = glob.glob(outFil)
        if (len(a) > 0):
            os.remove(a[0])
        fos=jio.FileOutputStream(outFil)
        ept.WriteSpectrumAsEMSA1_0.write(theAvg,fos,ept.WriteSpectrumAsEMSA1_0.Mode.COMPATIBLE)
        fos.close()
        if(debug):
            print(outFil)
    
def listCalibrations(det):
    """listCalibrations(det)
    list the calibrations for the specified detector.
    Example:
    import dtsa2.jmGen as jmg
    det = findDetector("FEI FIB620 EDAX-RTEM")
    jmg.listCalibrations(det)
    """
    print det
    tm=ju.TreeMap()
    for cal in dt2.Database.getCalibrations(det.getDetectorProperties()):
        tm.put(cal.getActiveDate(), cal)
    for me in tm.entrySet():
        print me.getValue()

def deleteDetector(det):
    """deleteDetector(det)
    Delete a detector. Use listDetectors() to identify the detector.
    By N. Ritchie at NIST. Use with care.    Any spectra associated with this detector
    will become lost souls.
    Example:
    import dtsa2.jmGen as jmg
    jmg.deleteDetector(d1)
    """
    dt2.Database.deleteDetector(dt2.Database.findDetector(det.getDetectorProperties()))
    
def dumpMaterials():
    """dumpMaterials()
    Dumping all the standards currently defined in the database to a Python script
    written to the user's home directory. The Python script can be used to export the
    standards to DTSA-II on another computer. By N. Ritchie at NIST. Note: I keep a
    annotated version in the productionScripts directory which can be added to from
    this output and is maintained under version control
    Example:
    import dtsa2.jmGen as jmg
    jmg.dumpMaterials() 
    """
    filename = jl.System.getProperty("user.home")+"/stds.py"
    # Get a set containing all standards
    stds=dt2.Database.findAllStandards()
    # Build up the repopulate helper script in the variable 'text' 
    text = "# -*- coding: utf-8 -*-\n\n"
    # Define a helper method 'defineStd'
    text += "def defineStd(elms,qty,name,density=None):\n"
    text += "    c=epq.Composition(map(element,elms),qty,name)\n"
    text += "    if density:\n"
    text += "        c=epq.Material(c,epq.ToSI.gPerCC(density))\n"
    text += "    Database.addStandard(c)\n\n"
    # For each standard in the database
    for mat in stds:
        elms = "(%s,)" % (",".join(["\"%s\"" % elm.toAbbrev() for elm in mat.getElementSet()]))
        qty = "(%s,)" % (",".join(["%f" % mat.weightFraction(elm,False) for elm in mat.getElementSet()]))
        text += "defineStd(%s,%s,\"" % (elms, qty)
        text += mat.toString()
        text += "\""
        if isinstance(mat,epq.Material):
            text += ",%f)\n" % epq.FromSI.gPerCC(mat.getDensity())
        else:
            text += ")\n"
    # Write the result to a text file.
    print text
    f = codecs.open(filename, 'w', encoding='utf-8')
    f.write(text)
    f.close()

def matchDet(spc, det):
    """matchDet(spc, det)
    Map the spectrum to match the channel length of the desired detector,
    copying the key parameters for quantitative analysis.
    Example:
    import dtsa2.jmGen as jmg
    det = findDetector("FEI FIB620 EDAX-RTEM")
    out = jmg.matchDet(spc, det)
    display(out)
    """
    dt2.display(spc)
    end = det.channelCount
    props = spc.getProperties()
    e0 = props.getNumericProperty(epq.SpectrumProperties.BeamEnergy)
    nm = props.getTextProperty(epq.SpectrumProperties.SpectrumDisplayName)
    ts = props.getTimestampProperty(epq.SpectrumProperties.AcquisitionTime)
    cw = spc.getChannelWidth()
    zo = spc.getZeroOffset()
    lt = spc.liveTime()
    pc = spc.probeCurrent()
    cr = epq.SpectrumUtils.slice(spc, 0, end)
    dt2.DataManager.removeSpectrum(spc)
    dt2.clear()
    sp = epq.SpectrumUtils.toSpectrum(cw, zo, end, cr)
    rm = dt2.wrap(sp)
    props = rm.getProperties()
    props.setDetector(det)
    props.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
    props.setTextProperty(epq.SpectrumProperties.SpectrumDisplayName, nm)
    props.setTimestampProperty(epq.SpectrumProperties.AcquisitionTime, ts)
    rm.setProbeCurrent(pc)
    rm.setLiveTime(lt)
    return rm
    
def tabulateDetCalibrations(det, outPath):
    """tabulateDetCalibrations(det, outPath)
    Tabulate the calibrations for the specified detector to
    the specified output file
    Example:
    import dtsa2.jmGen as jmg
    det = findDetector("FEI CM20UT EDAX-RTEM")
    out = jmg.tabulateDetCalibrations(det, './fei-cm20ut-det-cal.csv')
    """
    print(det.getName())
    dp = det.getDetectorProperties()
    cals = dt2.Database.getCalibrations(dp)
    iCnt = 0
    lDa = []
    lCw = []
    lZo = []
    lRe = []
    print(len(cals))
    for cal in cals:
        ad = cal.getActiveDate()
        lDa.append(ad.toString())
        cw = cal.getChannelWidth()
        lCw.append(cw)
        zo = cal.getZeroOffset()
        lZo.append(zo)
        cp = cal.getProperties()
        res = cp.getNumericProperty(cp.Resolution)
        lRe.append(res)
        lin = cp.getTextProperty(cp.ResolutionLine)
        iCnt += 1
    # write out results
    f = open(outPath, 'w')
    strLine = 'active.date, channel.width.eV, zero.offset.eV, resolution.eV\n'
    f.write(strLine)
    for i in range(iCnt):
        strLine = "%s" % lDa[i] + ","
        strLine = strLine + "%.5f" % lCw[i]     + ","
        strLine = strLine + "%.5f" % lZo[i]    + ","
        strLine = strLine + "%.2f" % lRe[i]    + "\n"
        f.write(strLine)
    f.close()
    
def elapsedTime(startTime):
    """elapsedTime(startTime)
    A convenience funtion to compute and print the elapsed time. Useful
    for long scripts we want to run in sequence. startTime is the starting
    time.
    Example
    import time
    import dtsa2.jmGen as jmg
    start = time.time()
    # other commands...
    jmg.elapsedTime(start)
    """
    end = time.time()
    delta = (end-startTime)/60
    msg = "Total computation time required %.3f min" % delta
    print msg
    if (delta > 60):
        delta = delta / 60
        msg = "...or %.2f hrs..." % delta

def intCuK(spc):
    """intCuK(spc)
    A convenience function to generate the background-subtracted peak
    integral of CuK and return in counts/per sec. This is to help to
    standardize beam measurement a la Oxford. Returns a list with the
    mean cps and the standard deviation"""
    spc = spc.applyLLT()
    lt    = spc.liveTime()
    # use the Ka and Kb peaks
    int = spc.peakIntegral(7661, 9331)
    piMuCps = int.doubleValue()/lt
    piSdCps = int.uncertainty()/lt
    ret = [piMuCps, piSdCps]
    return ret
    
def intCuL(spc):
    """intCuL(spc)
    A convenience function to generate the background-subtracted peak
    integral of CuL and return in counts/per sec. This is to help to
    standardize beam measurement a la Oxford. Returns a list with the
    mean cps and the standard deviation"""
    spc = spc.applyLLT()
    lt    = spc.liveTime()
    # use the Ka and Kb peaks
    int = spc.peakIntegral(633, 1147)
    piMuCps = int.doubleValue()/lt
    piSdCps = int.uncertainty()/lt
    ret = [piMuCps, piSdCps]
    return ret
