"""
test-mc3-simulate-pyrex.py
"""


import os, shutil
import dtsa2.mcSimulate3 as mc3

start = time.time()

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
outDir = homDir + "/Documents/git/dtsa2Scripts/utility/output/"
rptDir = wrkDir + '/test-mc3-simulate-K227 Results/'

DataManager.clearSpectrumList()

e0  = 20.0
det = findDetector("Probe")
mat = material("Fe2O3",5.0)

# Define the materials
k412 = epq.Material(epq.Composition([epq.Element.Ca,epq.Element.Fe,epq.Element.O, epq.Element.Al,epq.Element.Mg, epq.Element.Si],
                                    [0.1090, 0.0774, 0.4276, 0.0491, 0.1166, 0.2120] ),
                                    epq.ToSI.gPerCC(5.0))

spc = mc3.simulate(mat, det, e0, 120.0, True, 1000, True, True, {})
sName = "Bulk Fe2O3"
spc.rename(sName)
spc.display()

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
