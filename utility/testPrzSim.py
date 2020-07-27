"""
testPrzSim.py

"""

import os, shutil

start = time.time()

det = findDetector("Oxford p4 05eV 2K")
e0  = 15.     # keV
lt  = 100     # sec
pc  =   5.0   # nA
withPoisson = True

name = "SiO2"
sio2 = material("SiO2", density=2.65)
si  = material("Si", density=2.648)


"""
From the certificate for SRM 470

The values in the dictionary are mass fractions of the component
oxides.

name = 'K-412'
mix = mixture({ "SiO2"  : 0.4535,
                "MgO"   : 0.1933,
                "CaO"   : 0.1525,
                "FeO"   : 0.0996,
                "Al2O3" : 0.0927},
               density=2.600,
               name=name)
"""

# should not need to change below here

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
outDir = homDir + "/Documents/git/dtsa2Scripts/html"
rptDir = wrkDir + '/testPrzSim Results/'

dose = pc * lt  # na-sec"

DataManager.clearSpectrumList()



mix = material(sio2)
sp = epq.SpectrumProperties()
sp.setDetector(det)
sp.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
sp.setNumericProperty(epq.SpectrumProperties.FaradayBegin, pc)
sp.setNumericProperty(epq.SpectrumProperties.LiveTime, lt)
sp.setDetector(det)
spc = epq.SpectrumSimulator.Basic.generateSpectrum(sio2, sp, True)
print(dir(spc))

res = wrap(spc)
if withPoisson:
    res = wrap(epq.SpectrumUtils.addNoiseToSpectrum(res, 1.0))
res.rename("Simulated[%s,%0.1f keV,%0.2f nA\267s]" % (mix, e0, dose))
res.display()

ss = epq.SpectrumSimulator.Basic.shellSet(mix, sp)
print(ss)

ints = epq.SpectrumSimulator.Basic.computeIntensities(mix, sp, ss)
print(ints)

# print(ints["O K-L2"])
print(ints.keys())
print(ints.values())

intO = 0.0
intSi = 0.0

for k in ints.keys():
    v = ints[k]
    s = k.toString()
    print (s, type(s))
    x = s.split(" ")[0]
    if (x == 'O'):
        intO += v
    if (x == 'Si'):
        intSi += v

print(intO, intSi)

ratioOtoSi = round(intO/intSi, 5)
print(ratioOtoSi)


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
