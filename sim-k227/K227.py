# K227.py
s=dict(Database.getStandards())
mat = material("K227")
props = epq.SpectrumProperties()
props.addAll(d2.getProperties())
props.setNumericProperty(epq.SpectrumProperties.BeamEnergy,25.0)
xpp = epq.CorrectionAlgorithm.XPP
xrts=[transition("Si K-L3"),transition("O K-L3"),transition("Pb L3-M5"),transition("Pb M5-N7")]
mac=epq.MassAbsorptionCoefficient.Chantler2005
for xrt in xrts:
	matMac=mac.computeWithUncertaintyEstimate(mat,xrt)
	print "%s\t%f\t%f" % (xrt,matMac.doubleValue(),matMac.uncertainty())
	for name, idx in s.iteritems():
		std=material(name)
		elm = xrt.getElement()
		if std.containsElement(elm):
			stdMac=mac.computeWithUncertaintyEstimate(std,xrt)
			print "%s\t%f\t%f\t%f\t%f" % (std,  std.weightFraction(elm,False), xpp.u_chi(std,mat,xrt,props), stdMac.doubleValue(),stdMac.uncertainty())

"""
2021-06-15  J. R. Minter
----------

Running /Users/jrminter/Documents/git/dtsa2scripts/K227/K227.py
Si K-L3 136.550064   11.874185
O K-L3  985.465932  167.903408
Pb L3-M5  8.156100    0.783112
Pb M5-N7 86.380349    6.496248

Elapse: 0:00:00.4
"""

