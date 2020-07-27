# mc3SimulateBaTiO3.py
#
# Simulate a bulk barium titanante spectrum. This is as simple as it
# gets. 
#
# Date        Ver  Who  Notes
# 2015-05-21 0.90  JRM  Initial example. Verified with Iona v.2015-05-01
#
import dtsa2.mcSimulate3 as mc3

det         = findDetector("Probe") # DTSA-II default detector, use yours
e0          = 15   # keV
nTraj       = 1000 # electrons
dose        = 150  # nA*sec
withPoisson = True # Add Poisson noise
sf          = True # Secondary Florescence
bf          = True # Brehmstrahlung Florescence
 
mat = epq.Material(epq.Composition([epq.Element.Ba, epq.Element.Ti, epq.Element.O],
                                   [0.2058,         0.2053,         0.5889]      ),
                                    epq.ToSI.gPerCC(6.02))
mat.setName("BaTiO3")

spc = mc3.simulate(mat, det, e0, dose, withPoisson, nTraj) #  sf=False, bf=False, xtraParams={})

display(spc)
