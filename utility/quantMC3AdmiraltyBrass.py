"""
quantMC3AdmiraltyBrass.py

Do a MC3 simulatiion of Admiralty Brass and standards at 20 kV and try
a quant 

Elapse: 1:21:24.9 for 20,000 traj

<table class="leftalign">
<tr><th colspan="4">Admiralty Brass-at-20-kV 20,000 traj</th></tr>
	<tr><th>Element</th><th>Mass<br/>Fraction</th><th>Mass Fraction<br/>(normalized)</th><th>Atomic<br/>Fraction</th></tr>
	<tr><td>Cu</td><td>0.6899±0.0012</td><td>0.6924±0.0012</td><td>0.7017±0.0012</td></tr>
	<tr><td>Zn</td><td>0.2961±0.0010</td><td>0.2972±0.0010</td><td>0.2926±0.0010</td></tr>
	<tr><td>Sn</td><td>0.0104±0.0002</td><td>0.0104±0.0002</td><td>0.0057±0.0001</td></tr>
</table>

"""

import os, shutil, time
import dtsa2.dt2Conv as dt2c
import dtsa2.mcSimulate3 as mc3

start = time.time()


"""
From:
text and here:
https://www.thoughtco.com/common-brass-alloys-and-their-uses-603706

densities
https://www.angstromsciences.com/density-elements-chart

The values in the dictionary are mass fractions of the component
oxides.
"""

name = 'Admiralty Brass'
rho = 0.69*8.86 + 0.30*7.13 + 0.01*7.31
admb = mixture({ "Cu" : 0.69,
                 "Zn" : 0.30,
                 "Sn" : 0.01
                },
                density=rho,
                name=name)


cup = material("Cu2O", density=6.1)
cup.setName("Cuprite")

wil = material("Zn2SiO4", density=4.05)
wil.setName("Willemite")

cas = material("SnO2", density=6.9)
cas.setName("Cassiterite")

cu = material("Cu", density=8.96)
cu.setName("Copper")

zn = material("Zn", density=7.13)
zn.setName("Zinc")




det  = findDetector("Oxford p4 05eV 4K")
nTraj    = 20000     # trajectories
lt       =   100     # sec
pc       =     1.0   # nA
e0       =    20.0   # kV

dose = pc * lt  # na-sec"

# should not need to change below here

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
outDir = homDir + "/Documents/git/dtsa2Scripts/html"
rptDir = wrkDir + '/quantMC3AdmiraltyBrass Results/'
simDir = homDir + "/Documents/git/dtsa2Scripts/sim-Admiralty-Brass"

fiPath = outDir + "/%s-quant.html" % (name)
dt2c.ensureDir(simDir)

DataManager.clearSpectrumList()

def simSpc(mat, e0, det, dose, ntraj, simDir):
	"""
    simSpc(mat, e0, det, dose, ntraj, simDir)

    simulate a spectrum from a material and write it out.
	"""
	xrts = []
	trs = mc3.suggestTransitions(mat, e0)
	for tr in trs:
		xrts.append(tr)

	xtraParams={}
	xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))
	xtraParams.update(mc3.configureOutput(simDir))

	spc = mc3.simulate(mat, det, e0, dose, True,
                       nTraj, True, True, xtraParams)
	fmtS = "%s-at-%g-kV"
	sName = fmtS % (mat.getName(), e0)
	spc.rename(sName)
	spc.display()
	fi =  simDir + "/"
	fi += sName
	fi += "-%g-Traj.msa" % (nTraj)
	spc.save(fi)
	return spc

unkS = simSpc(admb, e0, det, dose, nTraj, simDir)

d1 = time.time()
delta = (d1 -start)/60
msg = "Unknown sim required %.3f min" % delta
print msg

cupS = simSpc(cup,  e0, det, dose, nTraj, simDir)

d2 = time.time()
delta = (d2 - d1)/60
msg = "Cuprite sim required %.3f min" % delta
print msg


wilS = simSpc(wil,  e0, det, dose, nTraj, simDir)
d3 = time.time()
delta = (d3 - d2)/60
msg = "Willemite sim required %.3f min" % delta
print msg

casS = simSpc(cas,  e0, det, dose, nTraj, simDir)
d4 = time.time()
delta = (d4 - d3)/60
msg = "Cassiterite sim required %.3f min" % delta
print msg


# cuS  = simSpc(cu,   e0, det, dose, nTraj, simDir)
# znS  = simSpc(zn,   e0, det, dose, nTraj, simDir)

ret = quantify(unkS, {"Cu" : cupS, "Zn" : wilS, "Sn" : casS },
               {},
               preferred=(),
               elmByDiff=None,
               oByStoic=False,
               oxidizer=None,
               extraKRatios=None,
               fiat={})

print(ret)
print(dir(ret))

print(ret.getUnknown())

a = ret.getComposition()
ht = a.toHTMLTable()

# f=open(fiPath, 'w')
print(ht)
# f.close()


print(ret.composition)

display(ret.getResidual())









# clean up cruft
shutil.rmtree(rptDir)
print "Done!"

end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg

