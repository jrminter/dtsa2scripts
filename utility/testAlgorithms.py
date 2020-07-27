# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

testAlgorithms.py

  Date      Who  Comment
----------  ---  -----------------------------------------------
2017-04-12  JRM  Figuring out mass-thickness (rho*z). Initial.


"""

import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQTools as ept
import gov.nist.microanalysis.Utility as epu
# This next two lines of weirdness eliminate a "no module named" error.  Why????
import sys as sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte", "MonteCarloSS", None)
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.microanalysis.dtsa2 as dt2
import gov.nist.microanalysis.EPQLibrary.Detector as epd
import java.lang as jl
import java.io as jio
import java.util as ju
import jarray
import javax.imageio as ii
import os, shutil

gitDir  = os.environ['GIT_HOME']
relPrj  = "/dtsa2Scripts/utility"
prjDir  = gitDir + relPrj
rptDir  = prjDir + '/testAlgorithms Results/'


e0 = 15 # kV

# Define the material k412
# Density from Probe Software, compo from NIST spectrum
k412 = epq.Material(epq.Composition([epq.Element.O,
                                     epq.Element.Mg,
                                     epq.Element.Al,
                                     epq.Element.Ca,
                                     epq.Element.Si,
                                     epq.Element.Fe],
                                    [ 0.427580,
                                      0.116567,
                                      0.049062,
                                      0.211982,
                                      0.108990,
                                      0.077420]
                                      ),
                                    epq.ToSI.gPerCC(2.600))
k412.setName("K412")

au = epq.Material(epq.Composition([epq.Element.Au], [1.0]), epq.ToSI.gPerCC(19.30))
au.setName("Au")

def ComputeElectronRange(mat, e0):
    """ComputeElectronRange(mat, e0)
    Compute the electron range for a material using two
    algorithms.

    Parameters
    ----------
    mat - epq.Material with a density
        the material to compute
    e0 - float
        the kV

    Returns
    -------
    none - writes to screen
    """

    # The result is in meters * (kg/meter^3)
    matName = mat.getName()
    rho = mat.getDensity()
    strOut = "Electron range computations for %s with density %g g/cm3" % (matName, 0.001*rho)
    print(strOut)
    rhoZmax = epq.ElectronRange.KanayaAndOkayama1972.compute(mat, epq.ToSI.keV(e0))
    zMax = 1.0e6*rhoZmax / rho
    rhoZmax *= 100.
    strOut = "Kanaya and Okayama (1972) Electron range: %.4f mg/cm^2 or %.4f um" % (rhoZmax, zMax)
    print(strOut)
    rhoZmax = epq.ElectronRange.Pouchou1991.compute(mat, epq.ToSI.keV(e0))
    zMax = 1.0e6*rhoZmax / rho
    rhoZmax *= 100.
    strOut = "Pouchou (1991)            Electron range: %.4f mg/cm^2 or %.4f um\n" % (rhoZmax, zMax)
    print(strOut)

ComputeElectronRange(k412, e0)
ComputeElectronRange(au, e0)

shutil.rmtree(rptDir)
print "Done!"