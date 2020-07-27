"""
testCalcKF.py

"""

import os, shutil
import dtsa2.dt2Conv as dt2c

start = time.time()

det = findDetector("Oxford p4 05eV 2K")
e0  = 15.     # keV
lt  = 100     # sec
pc  =   5.0   # nA
withPoisson = True

sp = epq.SpectrumProperties()
sio2 = material("SiO2", density=2.65)
si  = material("Si", density=2.648)
sio = material("SiO", density=2.13)

ver = False


# should not need to change below here

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
outDir = homDir + "/Documents/git/dtsa2Scripts/html"
rptDir = wrkDir + '/testCalcKF Results/'

dose = pc * lt  # na-sec"

DataManager.clearSpectrumList()

def calcEdsKF(comp, det, e0, alg=epq.XPP1991(),
              mac=epq.MassAbsorptionCoefficient.Chantler2005,
              xtra=epq.SpectrumProperties(),
              stds=None, verbose=True):
	"""
   calcEdsKF(comp, det, e0, alg=epq.XPP1991(),
             mac=epq.MassAbsorptionCoefficient.Chantler2005,
             xtra=epq.SpectrumProperties(),
             stds=None, verbose=True)

    Calculate EDS K-factors for a composition, otionally using
    standards. This can be helful for analyses using compound standards.

    Parameters
    ----------
    comp: epq.Composition
        The composition for K-factor computation
    det: epq.Library.EDSDetector
        The detector to specify key properties
    e0: float
        The beam energy in kV
    alg: algorithm for phirhoz (epq.XPP1991())
        epq.PAP1991() is another good choice
    mac: epq.MassAbsorptionCoefficient.Chantler2005
        Algorithm for mass absorption
    xtra: spectrum properties (epq.SpectrumProperties())
        Spectrum properties to pass detector
    stds: A dictionary of standards to use if desired (None)
        Standards. If you don't specify stds then pure elements are
        assumed. Otheerwise a dictionary mapping an element to a \
        composition Example: {"Si": sio2} where sio2 is a composition
    verbose: boolean (True)
    	print sll output

    Returns
    -------
    (lTrans, lK) - a list of transitions and a list of k-factors

    Example
    -------
    ver = False
    det = findDetector("Oxford p4 05eV 2K")
    e0  = 15.
    sp = epq.SpectrumProperties()
    sio2 = material("SiO2", density=2.65)
    si  = material("Si", density=2.648)
    resPAP = calcEdsKF(si, det, e0, epq.PAP1991(),
             epq.MassAbsorptionCoefficient.Chantler2005,
             xtra=sp, stds={"Si": sio}, verbose=ver)
    print(resPAP)
"""
	lK =[]
	lTrans = majorTransitionSets(det, comp, e0, 0.01)
	if isinstance(comp, epq.ISpectrumData):
		cp = comp.getProperties()
		t = cp.getCompositionWithDefault(epq.SpectrumProperties.MicroanalyticalComposition, None)
		if not t:
		  t = cp.getCompositionWithDefault(epq.SpectrumProperties.StandardComposition, None)
		comp = t
	comp = material(comp)
	if (comp <> None) and (comp.getElementCount() > 0):
		oldStrat = epq.AlgorithmUser.getGlobalStrategy()
		s = epq.Strategy()
		s.addAlgorithm(epq.MassAbsorptionCoefficient, mac)
		epq.AlgorithmUser.applyGlobalOverride(s)
		props = epq.SpectrumProperties()
		props.setDetector(det)
		props.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
		props.addAll(xtra)
		if verbose:
			print "Material\t%s" % (comp.descriptiveString(0))
			print "Detector\t%s" % det
			print "Algorithm\t%s" % alg.getName()
			print "MAC\t%s" % alg.getAlgorithm(epq.MassAbsorptionCoefficient).getName()
			print "E0\t%g keV" % props.getNumericProperty(epq.SpectrumProperties.BeamEnergy)
			print "Take-off\t%g%s" % (jl.Math.toDegrees(epq.SpectrumUtils.getTakeOffAngle(props)), epq.SpectrumProperties.TakeOffAngle.getUnits())
		for sp in xtra.getPropertySet():
			if verbose:
				print "%s\t%s" % (sp, xtra.getTextProperty(sp))
		if stds:
			conv = {}
			for z, c in stds.iteritems():
				conv[element(z)] = material(c)
			stds = conv
		if verbose:
			print "\n%-15s\tStandard\tEnergy\t ZAF\t  Z\t  A\t  F\tk-ratio" % "IUPAC"
		for xrts in majorTransitionSets(det, comp, e0, 0.01):
			z, a, f, zaf, w = 0.0, 0.0, 0.0, 0.0, 0.0
			elm = xrts.getElement()
			std = (stds.get(elm) if stds else None)
			ww = (std.weightFraction(elm, False) if std else 1.0)
			nComp = comp
			if comp.weightFraction(elm, False) < 1.0e-8:
				 nComp = epq.Composition(comp)
				 nComp.addElement(elm, 1.0e-8)
			for xrt in xrts:
				if epq.FromSI.keV(xrt.getEdgeEnergy()) > 0.9 * e0:
					 continue
				rzaf = (alg.relativeZAF(nComp, xrt, props, std) if std else alg.relativeZAF(nComp, xrt, props))
				wgt = xrt.getWeight(epq.XRayTransition.NormalizeFamily)
				z = z + wgt * rzaf[0] 
				a = a + wgt * rzaf[1] 
				f = f + wgt * rzaf[2] 
				zaf = zaf + wgt * rzaf[3]
				w = w + wgt
				eTr = epq.FromSI.keV(xrt.getEnergy())
			if w < 1.0e-10:
				continue
			z, a, f, zaf = z / w, a / w, f / w, zaf / w
			k = zaf * nComp.weightFraction(elm, False) / ww
			lK.append(k)
			if verbose:
				print "%-15s\t%s\t%2.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.6f" % (xrts, (std if std else "Pure %s" % elm.toAbbrev()), eTr, zaf, z, a, f, k)
	# Restore previous algorithm preferences
	epq.AlgorithmUser.applyGlobalOverride(oldStrat)
	return(lTrans, lK)

print("\n\n")
print("PAP1991-Chantler2005 SiO2-SiO")
print("\n")

resPAP = dt2c.calcEdsKF(sio2, det, e0, epq.PAP1991(),
             epq.MassAbsorptionCoefficient.Chantler2005,
             xtra=sp, stds={"O": sio, "Si": sio}, verbose=ver)
print(resPAP)

print("\n\n")
print("XPP1991-Chantler2005 Si-SiO")
print("\n")

resPAP = dt2c.calcEdsKF(si, det, e0, epq.PAP1991(),
             epq.MassAbsorptionCoefficient.Chantler2005,
             xtra=sp, stds={"Si": sio}, verbose=ver)
print(resPAP)

print("\n\n")
print("PAP1991-Chantler2005 Si-SiO2")
print("\n")

resPAP = dt2c.calcEdsKF(si, det, e0, epq.PAP1991(),
             epq.MassAbsorptionCoefficient.Chantler2005,
             xtra=sp, stds={"Si": sio2}, verbose=ver)
print(resPAP)


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
