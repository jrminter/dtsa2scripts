# -*- coding: utf-8 -*-
#        1         2         3         4         5         6         7 |
# 3456789012345678901234567890123456789012345678901234567890123456789012
#
# jmMC3.py
#
# Wrapper scripts for MC3 Monte Carlo simulations with DTSA-II
# Developed for Iona version 2015-07-01 
#
#    Modifications
#      Date   Who   Ver     What
# ----------  ---  ------  ---------------------------------------------
# 2015-07-13  JRM  0.0.10  Initial prototype: biLayerLineOnSubstrate, 
#                          triLayerLineOnSubstrate, and
#                          simulateBulkStandard.
# 2015-07-14  JRM  0.0.11  Added simLineInMatrix and formatted.
# 2015-07-23  JRM  0.0.12  Added simLineInMatrixLimScan
# 2015-10-02  JRM  0.0.13  Added coatedSphereOnFilm
# 2015-10-03  JRM  0.0.14  Added coatedOverBlock
# 2016-05-14  JRM  0.0.15  Updated for Jupiter
# 2017-06-18  JRM  0.0.16  Added simCarbonCoatedStd
# 2017-06-20  JRM  0.0.17  Added simBulkStd to sim uncoated standard
# 2017-06-22  JRM  0.0.18  Added simCtdOxOnSi
# 2018-10-14  JRM  0.0.19  Added getSpecPath

__revision__ = "$Id: jmMC3.py John Minter $"
__version__ = "0.0.19"

import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3, FluorescenceXRayGeneration3, XRayTransport3", None)
import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQLibrary.Detector as epd
import gov.nist.microanalysis.Utility as epu
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.microanalysis.NISTMonte.Gen3 as nm3
import gov.nist.microanalysis.EPQTools as et
import dtsa2.mcSimulate3 as mc3
import dtsa2.jmGen as jmg
import dtsa2 as dt2
import java.util as jutil
import java.io as jio
import java.nio.charset as cs
import os
import glob
import shutil

if 'defaultXtraParams' not in globals():
   defaultXtraParams = {}
if 'defaultBremFluor' not in globals():
   defaultBremFluor = False
if 'defaultCharFluor' not in globals():
   defaultCharFluor = False
if 'defaultNumTraj' not in globals():
   defaultNumTraj = 1000
if 'defaultDose' not in globals():
   defaultDose = 120.0

def getSpecPath(baseName, baseDir, e0, nTraj):
    """
    getSpecPath(baseName, baseDir, e0, nTraj)

    Generate a file path for Monte Carlo simulated spectrum

    Parameters
    ----------
    baseName - a string. The base name for the simulation,
               Example:  "bulkC"
    baseDir  - a string. The path to the directory to write the spectrum
               Example:  "C:/username/Documents/work/project"
                         Note: no path separator!
    e0       - a number. The voltage (kV) for the simulation
               Example:  15
    nTraj    - a number. The number of trajectories for the simulation
               Example:  20000

    Returns
    -------
    path     - A string. The path to the file to write

    Example
    -------
    import dtsa2.jmMC3 as jm3
    e0 = 15
    nTraj = 20000
    det = findDetector("Oxford p4 05eV 4K")
    c = material("C", density=2.266)
    a = jm3.simBulkStd(c, det, e0, nTraj, 100, 1.0, False)
    a.display()
    fi = jm3.getSpecPath("bulkC", "C:/username/Documents/work/project", e0, nTraj)
    a.save(fi)

    """
    sName = "%s-%g-kV-%g-Traj" % (baseName, e0, nTraj)
    sPath = "%s/%s.msa" % (baseDir, sName)
    return sPath

def simBulkStd(mat, det, e0, nTraj, lt=100, pc=1.0, ctd=True):
    """simBulkStd(mat, det, e0, nTraj, lt=100, pc=1.0)

    Use mc3 simulation to simulate an uncoated standard specimen

    Parameters
    ----------
    mat - a dtsa material.
        Note the material must have an associated density. It should have a useful name.
    det - a dtsa detector
        Here is where we get the detector properties and calibration
    e0 - float
        The accelerating voltage in kV
    nTraj - integer
        The number of trajectories to run
    lt - integer (100)
        The live time (sec)
    pc - float (1.0)
        The probe current in nA
    ctd - Boolean (True) - is C coated


    Returns
    -------
    sim - DTSA scriptable spectrum 
        The simulated standard spectrum
    
    Example
    -------
    import dtsa2.jmMC3 as jm3
    det = findDetector("Oxford p4 05eV 2K")
    cu = material("Cu", density=8.92)
    a = jm3.simBulkStd(cu, det, 20.0, 100, 100, 1.0)
    a.display()

    """
    dose = pc * lt  # na-sec"
    xrts = []

    trs = mc3.suggestTransitions(mat, e0)
    for tr in trs:
        xrts.append(tr)

    xtraParams={}
    xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))

    sim = mc3.simulate(mat, det, e0, dose, withPoisson=True, nTraj=nTraj,
                       sf=True, bf=True, xtraParams=xtraParams)

    sName = "%s-%g-kV" % (mat, e0)
    sim.rename(sName)
    sim.setAsStandard(mat)
    return sim

def simCtdOxOnSi(det, e0, nTraj, lt=100, pc=1.0, tox = 10.0, tc=20.0):
    """
    simCtdOxOnSi(det, e0, nTraj, lt=100, pc=1.0, tox = 10.0, tc=20.0)

    Use mc3 multilayer simulation to simulate a C-ctd silicon specimen
    with a native oxide layer.

    det - a dtsa detector
        Here is where we get the detector properties and calibration
    e0 - float
        The accelerating voltage in kV
    nTraj - integer
        The number of trajectories to run
    lt - integer (100)
        The live time (sec)
    pc - float (1.0)
        The probe current in nA
    tox - float (10.0)
        The thickness of the native oxide in nm
    tc - float (20.0)
        C thickness in nm


    Returns
    -------
    sim - DTSA scriptable spectrum 
        The simulated spectrum
    
    Example
    -------
    import dtsa2.jmMC3 as jm3
    det = findDetector("Oxford p4 05eV 2K")
    a = jm3.simCtdOxOnSi(det, 3.0, 100, 100, 1.0, 10.0, 20.0)
    a.display()

    """
    c    = dt2.material("C", density=2.1)
    si   = dt2.material("Si", density=2.329)
    sio2 = dt2.material("SiO2", density=2.65)

    dose = pc * lt  # na-sec"

    layers = [ [   c, tc*1.0e-9],
               [sio2, tox*1.0e-9],
               [si, 50.0e-6] ]
    xrts = []

    trs = mc3.suggestTransitions(c, e0)
    for tr in trs:
       xrts.append(tr)

    trs = mc3.suggestTransitions(sio2, e0)
    for tr in trs:
        xrts.append(tr)

    xtraParams={}
    xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))

    sim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                        dose=dose, sf=True, bf=True, xtraParams=xtraParams)
    sName = "%g-nm-C-on-%g-nm-SiO2-on-Si-%g-kV-%g-Traj" % (tc, tox, e0, nTraj)
    sim.rename(sName)
    return sim



def simCarbonCoatedStd(mat, det, e0, nTraj, lt=100, pc=1.0, tc=20.0):
    """simCarbonCoatedStd(mat, det, e0, nTraj, lt=100, pc=1.0, tc=20.0)

    Use mc3 multilayer simulation to simulate a C-ctd standard specimen

    Parameters
    ----------
    mat - a dtsa material.
        Note the material must have an associated density. It should have a useful name.
    det - a dtsa detector
        Here is where we get the detector properties and calibration
    e0 - float
        The accelerating voltage in kV
    nTraj - integer
        The number of trajectories to run
    lt - integer (100)
        The live time (sec)
    pc - float (1.0)
        The probe current in nA
    tc - float (20.0)
        C thickness in nm


    Returns
    -------
    sim - DTSA scriptable spectrum 
        The simulated standard spectrum
    
    Example
    -------
    import dtsa2.jmMC3 as jm3
    det = findDetector("Oxford p4 05eV 2K")
    mgo = material("MgO", density=3.58)
    a = jm3.simCarbonCoatedStd(mgo, det, 20.0, 100, 100, 1.0, 20.0)
    a.display()

    """
    dose = pc * lt  # na-sec"
    c = dt2.material("C", density=2.1)
    layers = [ [  c, tc*1.0e-9],
               [mat, 50.0e-6]
             ]
    xrts = []

    trs = mc3.suggestTransitions(c, e0)
    for tr in trs:
       xrts.append(tr)

    trs = mc3.suggestTransitions(mat, e0)
    for tr in trs:
        xrts.append(tr)

    xtraParams={}
    xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))

    sim = mc3.multiFilm(layers, det, e0, withPoisson=True, nTraj=nTraj,
                        dose=dose, sf=True, bf=True, xtraParams=xtraParams)
    sName = "%g-nm-C-on-%s-%g-kV-%g-Traj" % (tc, mat, e0, nTraj)
    sim.rename(sName)
    sim.setAsStandard(mat)
    return sim




def coatedOverBlock(mat, height, width, coating, thickness, substrate, det, e0=20.0, withPoisson=True, nTraj=defaultNumTraj, dose=defaultDose, sf=defaultCharFluor, bf=defaultBremFluor, xtraParams=defaultXtraParams):
    """coatedOverBlock(mat, height, width, coating, thickness, substrate, det, e0=20.0, withPoisson=True, nTraj=defaultNumTraj, dose=defaultDose, sf=defaultCharFluor, bf=defaultBremFluor, substrate=None, xtraParams={})
    Monte Carlo simulate a spectrum from a block shaped particle of the specified material (mat) and height (z in m) and width (x and y in m). \
    The block and subtrate is coated in a material 'coating' of the specified thickness which fully encapsulates the particle and covers the substrate too."""
    def buildBlock(monte, chamber, origin, buildParams):
        height = buildParams["Height"]
        width = buildParams["Width"]
        subMat = buildParams["Substrate"]
        mat = buildParams["Material"]
        coating = buildParams["Coating"]
        thickness = buildParams["Thickness"]
        
        coatedCube = nm.MultiPlaneShape.createBlock([width + 2.0 * thickness, width + 2.0 * thickness, height + thickness], epu.Math2.plus(origin, [0.0, 0.0, 0.5*(height + thickness)]), 0.0, 0.0, 0.0)
        sr1 = monte.addSubRegion(chamber, coating, coatedCube)
        cube = nm.MultiPlaneShape.createBlock([width, width, height], epu.Math2.plus(origin, [0.0, 0.0, thickness+0.5*height]), 0.0, 0.0, 0.0)
        monte.addSubRegion(sr1, mat, cube)
        sideSlabWidth = 2.5*1.0e-6 - (thickness + 0.5*width)
        sideSlabDims = [sideSlabWidth, 5.*1.0e-6, thickness]
        leftSidePt = epu.Math2.plus(origin, [0.5*(width+sideSlabWidth), 0.0, thickness+height])
        leftSide = nm.MultiPlaneShape.createBlock(sideSlabDims, leftSidePt, 0.0, 0.0, 0.0)
        monte.addSubRegion(chamber, coating, leftSide)
        rightSidePt = epu.Math2.plus(origin, [-0.5*(width+sideSlabWidth), 0.0, thickness+height])
        rightSide = nm.MultiPlaneShape.createBlock(sideSlabDims, rightSidePt, 0.0, 0.0, 0.0)
        monte.addSubRegion(chamber, coating, rightSide)
        fbSlabDims  = [width, sideSlabWidth, thickness]
        frontSidePt = epu.Math2.plus(origin, [0.0, 0.5*(width+sideSlabWidth), thickness+height])
        frontSide = nm.MultiPlaneShape.createBlock(fbSlabDims, frontSidePt, 0.0, 0.0, 0.0)
        monte.addSubRegion(chamber, coating, frontSide)
        backSidePt = epu.Math2.plus(origin, [0.0, -0.5*(width+sideSlabWidth), thickness+height])
        backSide = nm.MultiPlaneShape.createBlock(fbSlabDims, backSidePt, 0.0, 0.0, 0.0)
        monte.addSubRegion(chamber, coating, backSide)
        # no substrate - don't want film under block...
        # monte.addSubRegion(chamber, coating, nm.MultiPlaneShape.createFilm([0.0, 0.0, -1.0], epu.Math2.plus(origin, [0.0, 0.0, height + thickness]), thickness))
        monte.addSubRegion(chamber, subMat, nm.MultiPlaneShape.createSubstrate([0.0, 0.0, -1.0], epu.Math2.plus(origin, [0.0, 0.0, height + 2*thickness])))

    s1 = u"MC simulation of a [%0.2f,%0.2f,%0.2f] micron block of %s%s" % (width * 1.0e6, width * 1.0e6, height * 1.0e6, mat, (" on %s" % substrate if substrate else ""))
    s2 = u" coated with %0.2f microns of %s at %0.1f keV%s%s" % (thickness* 1.0e6, coating, e0, (" + CSF" if sf else ""), (" + BSF" if bf else ""))
    tmp = s1 + s2

    # tmp = u"MC simulation of a [%0.2f,%0.2f,%0.2f] micron block of %s%s coated with %0.2f microns of %s at %0.1f keV%s%s" % (width * 1.0e6, width * 1.0e6, height * 1.0e6, mat, (" on %s" % substrate if substrate else ""), coating, e0, (" + CSF" if sf else ""), (" + BSF" if bf else ""))
    params = {"Substrate": substrate, "Width" : width, "Height" : height, "Material" : mat, "Coating" : coating, "Thickness" : thickness}
    return mc3.base(det, e0, withPoisson, nTraj, dose, sf, bf, tmp, buildBlock, params, xtraParams)



def coatedSphereOnFilm(mat, radius, coating, cThick, film, fThick, det, e0=20.0, withPoisson=True, nTraj=100, dose=100, sf=True, bf=True,  xtraParams={}):
    """coatedSphereOnFilm(mat, radius, coating, cThick, film, fThick, det, e0=20.0, withPoisson=True, nTraj=100, dose=100, sf=True, bf=True,  xtraParams={})
Monte Carlo simulate a spectrum from a spherical particle of the specified material (mat) and radius (in m). \
on a film of material film and thickness fThick immediately \
below the particle."""
    if radius < 0.0:
         raise "The sphere radius must be larger than zero."
    if cThick < 0.0:
         raise "The coating thickness must be larger than zero."
    if fThick < 0.0:
         raise "The film thickness must be larger than zero."
    def buildSphere(monte, chamber, origin, buildParams):
        mat = buildParams["Material"]
        radius = buildParams["Radius"]
        coating = buildParams["Coating"]
        cThick = buildParams["Coating Thickness"]
        film = buildParams["Film"]
        fThick = buildParams["Film Thickness"]
        coatSphere = nm.Sphere(epu.Math2.plus(origin, [0.0, 0.0, radius + cThick]), radius + cThick)
        srC = monte.addSubRegion(chamber, coating, coatSphere)
        sphere = nm.Sphere(epu.Math2.plus(origin, [0.0, 0.0, radius + cThick]), radius)
        monte.addSubRegion(srC, mat, sphere)
        monte.addSubRegion(chamber, film, nm.MultiPlaneShape.createFilm([0.0, 0.0, -1.0], epu.Math2.plus(origin, [0.0, 0.0, 2.0 * radius]), fThick))
    tmp = u"MC simulation of a %0.3f micron sphere of %s coated with %0.3f microns of %s on %0.3f microns of %s at %0.1f keV%s%s" % (radius * 1.0e6, mat, cThick * 1.0e6, coating, fThick* 1.0e6, film, e0, (" + CSF" if sf else ""), (" + BSF" if bf else ""))
    return mc3.base(det, e0, withPoisson, nTraj, dose, sf, bf, tmp, buildSphere, {"Film": film, "Radius" : radius, "Material" : mat, "Coating" : coating, "Coating Thickness" : cThick, "Film Thickness" : fThick}, xtraParams)



def biLayerLineOnSubstrate(matTopL, matBotL, matSub, htTopL, htBotL, width, length, det, title, e0=20.0, withPoisson=True, nTraj=100, dose=120.0, sf=True, bf=True, xtraParams={}, bVerbose=False):
    """biLayerLineOnSubstrate(matTopL, matBotL, matSub, htTopL, htBotL,
    width, length, det, title, e0=20.0, withPoisson=True,
    nTraj=100, dose=120.0, sf=True,
    bf=True, xtraParams={}, bVerbose=False)

    Monte Carlo simulate a spectrum from a bilayer line on a substrate of
    the specified materials (matTopL, matBotL, matSub) and layer heights
    ( htTopL, htBotL; z in m ) and width and length (x and y in m)."""

    def buildBlock(monte, chamber, origin, buildParams):
        matSub = buildParams["Mat Subs"]
        matTopL = buildParams["Mat Top"]
        matBotL = buildParams["Mat Bot"]
        htTopL = buildParams["Hei Top"]
        htBotL = buildParams["Hei Bot"]
        width = buildParams["Width"]
        length = buildParams["Length"]
        monte.addSubRegion(chamber, matTopL, nm.MultiPlaneShape.createBlock([width, length, htTopL], epu.Math2.plus(origin, [0.0, 0.0, 0.5*htTopL]), 0.0, 0.0, 0.0))
        monte.addSubRegion(chamber, matBotL, nm.MultiPlaneShape.createBlock([width, length, htBotL], epu.Math2.plus(origin, [0.0, 0.0, 0.5*htBotL+htTopL]), 0.0, 0.0, 0.0))
        monte.addSubRegion(chamber, matSub, nm.MultiPlaneShape.createSubstrate([0.0, 0.0, -1.0], epu.Math2.plus(origin, [0.0, 0.0, htTopL+htBotL])))

    tmp = u"biLayerLineOnSub-%s" % title
    params = {"Mat Top" : matTopL, "Mat Bot" : matBotL, "Mat Subs": matSub, "Hei Top" : htTopL, "Hei Bot" : htBotL, "Width" : width, "Length" : length}
    if (bVerbose==120.0):
        print(params)
    return mc3.base(det, e0, withPoisson, nTraj, dose, sf, bf, tmp, buildBlock, params, xtraParams)


def triLayerLineOnSubstrate(matTopL, matMidL, matBotL, matSub, htTopL, htMidL, htBotL, width, length, det, title, e0=20.0, withPoisson=True, nTraj=100, dose=120.0, sf=True, bf=True, xtraParams={}, bVerbose=False):
    """triLayerLineOnSubstrate(matTopL, matMidL, matBotL, matSub, htTopL,
    htMidL, htBotL, width, length, det, title, e0=20.0, withPoisson=True,
    nTraj=100, dose=120.0, sf=True,
    bf=True, xtraParams={})

    Monte Carlo simulate a spectrum from a bilayer line on a substrate of
    the specified materials (matTopL, matMid, matBotL, matSub) and layer
    heights ( htTopL, htMidL,    htBotL; z in m ) and width and length
    (x and y in m)."""
    def buildBlock(monte, chamber, origin, buildParams):
        matSub = buildParams["Mat Subs"]
        matTopL = buildParams["Mat Top"]
        matMidL = buildParams["Mat Mid"]
        matBotL = buildParams["Mat Bot"]
        htTopL = buildParams["Hei Top"]
        htMidL = buildParams["Hei Mid"]
        htBotL = buildParams["Hei Bot"]
        width = buildParams["Width"]
        length = buildParams["Length"]
        monte.addSubRegion(chamber, matTopL, nm.MultiPlaneShape.createBlock([width, length, htTopL], epu.Math2.plus(origin, [0.0, 0.0, 0.5*htTopL]), 0.0, 0.0, 0.0))
        monte.addSubRegion(chamber, matMidL, nm.MultiPlaneShape.createBlock([width, length, htMidL], epu.Math2.plus(origin, [0.0, 0.0, 0.5*htMidL+htTopL]), 0.0, 0.0, 0.0))
        monte.addSubRegion(chamber, matBotL, nm.MultiPlaneShape.createBlock([width, length, htBotL], epu.Math2.plus(origin, [0.0, 0.0, 0.5*htBotL+htTopL+htMidL]), 0.0, 0.0, 0.0))
        monte.addSubRegion(chamber, matSub, nm.MultiPlaneShape.createSubstrate([0.0, 0.0, -1.0], epu.Math2.plus(origin, [0.0, 0.0, htTopL+htMidL+htBotL])))
    tmp = u"triLayerLineOnSub-%s" % title
    params = {"Mat Top" : matTopL, "Mat Mid" : matMidL, "Mat Bot" : matBotL, "Mat Subs": matSub, "Hei Top" : htTopL, "Hei Mid" : htMidL, "Hei Bot" : htBotL, "Width" : width, "Length" : length}
    if (bVerbose==120.0):
        print(params)
    return mc3.base(det, e0, withPoisson, nTraj, dose, sf, bf, tmp, buildBlock, params, xtraParams)

def simulateBulkStandard(mat, name, det, e0, lt, pc, withPoisson=True, nTraj=100, sf=True, bf=True, xtraParams={}):
    """simulateBulkStandard(mat, name, det, e0, lt, pc, withPoisson=True,
    nTraj=100, sf=True, bf=True, xtraParams={})"""
    std = mc3.simulate(mat, det, e0, lt*pc, withPoisson=True, nTraj=nTraj, sf=True, bf=True, xtraParams={})
    props=std.getProperties()
    props.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
    props.setNumericProperty(epq.SpectrumProperties.FaradayBegin, pc)
    props.setNumericProperty(epq.SpectrumProperties.FaradayEnd, pc)
    props.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
    std.setAsStandard(mat)
    return(std)

def simLineInMatrix(lin, linMat, blk, blkMat, nmLinWid, umBlock, nPts, trs, outDir, hdr, det, e0, lt, pc, withPoisson=True, nTraj=100, sf=True, bf=True, iDigits=5, bVerbose=False, xtraParams={}):
    """simLineInMatrix(lin, linMat, blk, blkMat, nmLinWid, umBlock, nPts,
    trs, outDir, hdr, det, e0, lt, pc, withPoisson=True, nTraj=nTraj,
    sf=True, bf=True, iDigits=5, bVerbose=False, xtraParams={})
    Simulate a line of width `nmLinWid' nm at the center of a block of
    `umBlock' microns. The line is of material `lin' with a name `linMat'.
    The block is of material `blk' with a name `blkMat'.
    We analyze an list `trs' of transitions, writing the K-ratios to a
    .csv file with a header `hdr'. We use the detector `det', voltage `e0'
    (kV) and live time `lt' sec and probe current  `pc' nA. This will
    compute the standard spectra, compute the spectra for nPts+1 from
    -nPts/2 ... 0 ...nPts/2 times the block size. It will then compute the
    K-ratios for each spectrum and write them to a file `name' in outDir
    with a header `hdr' that matches the transition order.
    """
    # order is order of trs..
    sc = 1.0e-6 # scale from microns to meters for positions
    lX = [] # an array for postions
    lKlin = [] # an array for the K-ratio of the line
    lKblk = [] # an array for the K-ratio of the block. Title correspond to hdr string
    # start clean
    dt2.DataManager.clearSpectrumList()
    # create the standards
    linStd = simulateBulkStandard(lin, linMat, det, e0, lt, pc, withPoisson=True, nTraj=nTraj, sf=True, bf=True, xtraParams={})
    dt2.display(linStd)
    blkStd = simulateBulkStandard(blk, blkMat, det, e0, lt, pc, withPoisson=True, nTraj=nTraj, sf=True, bf=True, xtraParams={})
    dt2.display(blkStd)
    lStd = {"El":dt2.element(linMat), "Spc":linStd}
    bStd = {"El":dt2.element(blkMat), "Spc":blkStd}
    stds = [lStd, bStd] # note: put the transitions in this order
    iCount = 0
    for x in range(-nPts/2, (nPts/2)+1, 1):
        xv = sc*x*umBlock/nPts
        lX.append(round(x*umBlock/nPts, iDigits))
        monte=nm.MonteCarloSS()
        monte.setBeamEnergy(epq.ToSI.keV(e0))
        # use a 1 nm probe
        beam=nm.GaussianBeam(1.0e-9) 
        monte.setElectronGun(beam)
        beam.setCenter([xv, 0.0,-0.05])
        # createBlock(double[] dims, double[] point, double phi, double theta, double psi)
        # createBlock - Create a block of:
        #          dimensions specified in dims,
        #          centered at point,
        #          then rotated by the euler angles phi, theta, psi.
        block = nm.MultiPlaneShape.createBlock([umBlock*1.0e-6, umBlock*1.0e-6, umBlock*1.0e-6],[0.0,0.0, 0.5*umBlock*1.0e-6],0.0,0.0,0.0)
        matrix = monte.addSubRegion(monte.getChamber(), blk, block)
        monte.addSubRegion(matrix, lin, nm.MultiPlaneShape.createBlock([1.0e-9*nmLinWid, umBlock*1.0e-6, umBlock*1.0e-6],[0.0, 0.0, 0.5*umBlock*1.0e-6],0.0,0.0,0.0))
        det.reset()
        # Add event listeners to model characteristic radiation
        chXR = nm3.CharacteristicXRayGeneration3.create(monte)
        xrel = nm3.XRayTransport3.create(monte, det, chXR)
        brXR = nm3.BremsstrahlungXRayGeneration3.create(monte)
        brem = nm3.XRayTransport3.create(monte, det, brXR)
        fxg3 = nm3.FluorescenceXRayGeneration3.create(monte, chXR)
        chSF = nm3.XRayTransport3.create(monte, det, fxg3)
        brSF = nm3.XRayTransport3.create(monte, det, nm3.FluorescenceXRayGeneration3.create(monte, brXR))
        # here is where we run the simulation
        monte.runMultipleTrajectories(nTraj)
        spec = det.getSpectrum((lt*pc*1.0e-9) / (nTraj * epq.PhysicalConstants.ElectronCharge))
        props = spec.getProperties()
        props.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
        props.setNumericProperty(epq.SpectrumProperties.FaradayBegin, pc)
        props.setNumericProperty(epq.SpectrumProperties.FaradayEnd, pc)
        props.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
        spcName = "x = %.3f um" % x
        epq.SpectrumUtils.rename(spec, spcName)
        spec = epq.SpectrumUtils.addNoiseToSpectrum(spec, 1.0)
        # display(spec)
        a = jmg.compKRs(spec, stds, trs, det, e0)
        iCount += 1
        print(iCount)
        lKlin.append(round(a[0], iDigits))
        lKblk.append(round(a[1], iDigits))

    basFile ="%gnm-%s-in-%gum-%s-%gkV-%g-Traj.csv" % (nmLinWid, linMat, umBlock, blkMat, e0, nTraj)
    strOutFile = outDir + "/" + basFile
    f=open(strOutFile, 'w')
    strLine = hdr + '\n'
    f.write(strLine)
    for i in range(iCount):
        strLine = "%.5f" % lX[i] + ","
        strLine = strLine + "%.5f" % lKlin[i] + ","
        strLine = strLine + "%.5f" % lKblk[i] + "\n" 
        f.write(strLine)  
    f.close()


def simLineInMatrix3(lin, linMat, blk, blkMat, nmLinWid, umBlock, nPts, trs, outDir, hdr, det, e0, lt, pc, withPoisson=True, nTraj=100, sf=True, bf=True, iDigits=5, bVerbose=False, xtraParams={}):
    """simLineInMatrix3(lin, linMat, blk, blkMat, nmLinWid, umBlock, nPts,
    trs, outDir, hdr, det, e0, lt, pc, withPoisson=True, nTraj=nTraj,
    sf=True, bf=True, iDigits=5, bVerbose=False, xtraParams={})
    Simulate a line of width `nmLinWid' nm at the center of a block of
    `umBlock' microns. The line is of material `lin' with a name `linMat'.
    The block is of material `blk' with a name `blkMat'.
    We analyze an list `trs' of transitions, writing the K-ratios to a
    .csv file with a header `hdr'. We use the detector `det', voltage `e0'
    (kV) and live time `lt' sec and probe current  `pc' nA. This will
    compute the standard spectra, compute the spectra for nPts+1 from
    -nPts/2 ... 0 ...nPts/2 times the block size. It will then compute the
    K-ratios for each spectrum and write them to a file `name' in outDir
    with a header `hdr' that matches the transition order.
    """
    # order is order of trs..
    sc = 1.0e-6 # scale from microns to meters for positions
    dose = lt*pc
    lX = [] # an array for postions
    lKlin = [] # an array for the K-ratio of the line
    lKblk = [] # an array for the K-ratio of the block. Title correspond to hdr string
    umLine = nmLinWid * 1.0e-3
    # start clean
    dt2.DataManager.clearSpectrumList()
    # create the standards
    linStd = simulateBulkStandard(lin, linMat, det, e0, lt, pc, withPoisson=withPoisson, nTraj=nTraj, sf=sf, bf=bf, xtraParams={})
    dt2.display(linStd)
    blkStd = simulateBulkStandard(blk, blkMat, det, e0, lt, pc, withPoisson=withPoisson, nTraj=nTraj, sf=sf, bf=sf, xtraParams={})
    dt2.display(blkStd)
    lStd = {"El":dt2.element(linMat), "Spc":linStd}
    bStd = {"El":dt2.element(blkMat), "Spc":blkStd}
    stds = [lStd, bStd] # note: put the transitions in this order
    iCount = 0
    for x in range(-nPts/2, (nPts/2)+1, 1):
        xv = sc*x*umBlock/nPts
        lX.append(round(x*umBlock/nPts, iDigits))
        xtraParams={}
        xtraParams.update(mc3.configureXRayAccumulators(trs, charAccum=sf, charFluorAccum=sf, bremFluorAccum=bf))
        xtraParams.update(mc3.configureOutput(outDir))
        xtraParams.update(mc3.configureBeam(xv, 0, -0.099, 1.0))
        spec = mc3.embeddedRectangle(lin, [umLine*sc, umBlock*sc, umBlock*sc], blk, 0, det, e0, withPoisson=withPoisson, nTraj=nTraj, dose=dose, sf=sf, bf=bf, xtraParams=xtraParams)
        props = spec.getProperties()
        props.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
        props.setNumericProperty(epq.SpectrumProperties.FaradayBegin, pc)
        props.setNumericProperty(epq.SpectrumProperties.FaradayEnd, pc)
        props.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
        spcName = "x = %.3f um" % x
        epq.SpectrumUtils.rename(spec, spcName)
        spec = epq.SpectrumUtils.addNoiseToSpectrum(spec, 1.0)
        # display(spec)
        a = jmg.compKRs(spec, stds, trs, det, e0)
        iCount += 1
        print(iCount)
        lKlin.append(round(a[0], iDigits))
        lKblk.append(round(a[1], iDigits))

    basFile ="%gnm-%s-in-%gum-%s-%gkV-%g-Traj.csv" % (nmLinWid, linMat, umBlock, blkMat, e0, nTraj)
    strOutFile = outDir + "/" + basFile
    f=open(strOutFile, 'w')
    strLine = hdr + '\n'
    f.write(strLine)
    for i in range(iCount):
        strLine = "%.5f" % lX[i] + ","
        strLine = strLine + "%.5f" % lKlin[i] + ","
        strLine = strLine + "%.5f" % lKblk[i] + "\n" 
        f.write(strLine)  
    f.close()

def simLineInMatrixLimScan(lin, linMat, blk, blkMat, nmLinWid, umBlock,nmScan, nPts, trs, outDir, hdr, det, e0, lt, pc, withPoisson=True, nTraj=100, sf=True, bf=True, iDigits=5, bVerbose=False, xtraParams={}):
    """simLineInMatrixLimScan(lin, linMat, blk, blkMat, nmLinWid, umBlock,
    nmScan, nPts, trs, outDir, hdr, det, e0, lt, pc, withPoisson=True,
    nTraj=nTraj, sf=True, bf=True, iDigits=5, bVerbose=False,
    xtraParams={})

    Simulate a line of width `nmLinWid' nm at the center of a block of
    `umBlock' microns. The line is of material `lin' with a name `linMat'.
    The block is of material `blk' with a name `blkMat'. We step a total
    distance of nmScan across the center of the line.

    We analyze an list `trs' of transitions, writing the K-ratios to a
    .csv file with a header `hdr'. We use the detector `det', voltage `e0'
    (kV) and live time `lt' sec and probe current  `pc' nA. This will
    compute the standard spectra, compute the spectra the scanned
    region. It will then compute the K-ratios for each spectrum and write
    them to a file `name' in outDir with a header `hdr' that matches the
    transition order.
    """
    # order is order of trs..
    sc = 1.0e-6 # scale from microns to meters for positions
    dose = lt*pc
    lX = [] # an array for postions
    lKlin = [] # an array for the K-ratio of the line
    lKblk = [] # an array for the K-ratio of the block. Title correspond to hdr string
    umLine = nmLinWid * 1.0e-3
    # start clean
    dt2.DataManager.clearSpectrumList()
    # create the standards
    linStd = simulateBulkStandard(lin, linMat, det, e0, lt, pc, withPoisson=withPoisson, nTraj=nTraj, sf=sf, bf=bf, xtraParams={})
    dt2.display(linStd)
    blkStd = simulateBulkStandard(blk, blkMat, det, e0, lt, pc, withPoisson=withPoisson, nTraj=nTraj, sf=sf, bf=sf, xtraParams={})
    dt2.display(blkStd)
    lStd = {"El":dt2.element(linMat), "Spc":linStd}
    bStd = {"El":dt2.element(blkMat), "Spc":blkStd}
    stds = [lStd, bStd] # note: put the transitions in this order
    iCount = 0
    for x in range(-nPts/2, (nPts/2)+1, 1):
        xPosNm = x * nmScan / nPts
        lX.append(round(xPosNm, iDigits))
        xtraParams={}
        xtraParams.update(mc3.configureXRayAccumulators(trs, charAccum=sf, charFluorAccum=sf, bremFluorAccum=bf))
        xtraParams.update(mc3.configureOutput(outDir))
        xtraParams.update(mc3.configureBeam(xPosNm*1.0e-09, 0, -0.099, 1.0))
        spec = mc3.embeddedRectangle(lin, [umLine*sc, umBlock*sc, umBlock*sc], blk, 0, det, e0, withPoisson=withPoisson, nTraj=nTraj, dose=dose, sf=sf, bf=bf, xtraParams=xtraParams)
        props = spec.getProperties()
        props.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
        props.setNumericProperty(epq.SpectrumProperties.FaradayBegin, pc)
        props.setNumericProperty(epq.SpectrumProperties.FaradayEnd, pc)
        props.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
        spcName = "x = %.3f um" % x
        epq.SpectrumUtils.rename(spec, spcName)
        spec = epq.SpectrumUtils.addNoiseToSpectrum(spec, 1.0)
        # display(spec)
        a = jmg.compKRs(spec, stds, trs, det, e0)
        iCount += 1
        print(iCount, xPosNm)
        lKlin.append(round(a[0], iDigits))
        lKblk.append(round(a[1], iDigits))

    basFile ="%gnm-%s-in-%gum-%s-%gkV-%g-Traj.csv" % (nmLinWid, linMat, umBlock, blkMat, e0, nTraj)
    strOutFile = outDir + "/" + basFile
    f=open(strOutFile, 'w')
    strLine = hdr + '\n'
    f.write(strLine)
    for i in range(iCount):
        strLine = "%.3f" % lX[i] + ","
        strLine = strLine + "%.5f" % lKlin[i] + ","
        strLine = strLine + "%.5f" % lKblk[i] + "\n" 
        f.write(strLine)  
    f.close()

def lineInMatrix(lin, blk, nmLinWid, umBlock, det, e0=20.0, withPoisson=True, nTraj=100, dose=120.0, sf=True, bf=True, xtraParams={}):
    """lineInMatrix(lin, blk, nmLinWid, umBlock, det,
    e0=20.0, withPoisson=True, nTraj=100, dose=120.0,
    sf=True, bf=True,
    xtraParams={}"""
    def buildBlock(monte, chamber, origin, buildParams):
        lin = buildParams["Line"]
        blk = buildParams["Block"]
        nmLinWid = buildParams["Width"]
        umBlock = buildParams["Size"]
        sc = 1.0e-6 # scale from microns to meters for positions
        # createBlock(double[] dims, double[] point, double phi, double theta, double psi)
        # createBlock - Create a block of:
        #          dimensions specified in dims,
        #          centered at point,
        #          then rotated by the euler angles phi, theta, psi.
        block = nm.MultiPlaneShape.createBlock([umBlock*sc, umBlock*sc, umBlock*sc],[0.0,0.0, 0.5*umBlock*sc],0.0,0.0,0.0)
        matrix = monte.addSubRegion(monte.getChamber(), blk, block)
        monte.addSubRegion(matrix, lin, nm.MultiPlaneShape.createBlock([1.0e-9*nmLinWid, umBlock*sc, umBlock*sc],[0.0, 0.0, 0.5*umBlock*sc],0.0,0.0,0.0))
    tmp = u"MC3-sim-%g-nm-%s-line-in-%g-um-%s-block-%0.1f-kV" % (nmLinWid, lin, umBlock, blk, e0)
    params = { "Line": lin, "Width" : nmLinWid, "Block" : blk, "Size" : umBlock }
    return (mc3.base(det, e0, withPoisson, nTraj, dose, sf, bf, tmp, buildBlock, params, xtraParams))

