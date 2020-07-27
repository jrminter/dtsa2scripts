# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

k227-quant.py

A reproducible example to simulate and quantify K227 at 5 kV

  Date      Who   Ver    Comment
----------  ---  -----   ------------------------------------------
2020-05-31  JRM  0.0.1   Initial test. There is an issue with 
                         Note: This fails...
2020-06-01  JRM          Note output only reports Pb.

Running /Users/jrminter/Documents/git/dtsa2Scripts/sim-k227-quant/sim-k227-quant.py
[O K-L3, Si K-L3, Si K-M3, Pb M5-N7, Pb M4-N6]
MC simulation of bulk k227 at 5.0 keV + CSF + BSF
/Users/jrminter/dat/k227-quant/spc_k227.msa
[O K-L3, Si K-L3, Si K-M3]
MC simulation of bulk SiO at 5.0 keV + CSF + BSF
/Users/jrminter/dat/k227-quant/spc_sio.msa
[Pb M5-N7, Pb M4-N6]
MC simulation of bulk Pb at 5.0 keV + CSF + BSF
/Users/jrminter/dat/k227-quant/spc_pb.msa
k227-5kV = [Pb(0.5684±0.0027 mass frac),Σ=0.5684±0.0027]
335.904999971
Elapse: 0:05:35.9

From script

k227-5kV = [O(0.0466±0.0002 mass frac),Pb(0.6029±0.0028 mass frac),Σ=0.6494±0.0029]

Try again with:
print(tmp.descriptiveString(True)) - No difference... Set back to False

Try oByStoic True, still just get Pb

Jython 2.7.1 (default:0df7adb1b397, Jun 30 2017, 19:02:43) 
[OpenJDK 64-Bit Server VM (Amazon.com Inc.)] on java11.0.4

My Prefs:
Correction Algorithm: PAP Full phi-rho-z
Mass Abs Coef NIST-Chandler 2005
Brhem and. disn: Acosta 2002 linear
Ioniz Cross Sec: Bote/Salvat 2008


If I run a quant manually from the GUI I get

Spectrum  Sum           O             Si            Pb
spc_k227  0.9656±0.0142 0.1619±0.0136 0.0941±0.0007 0.7096±0.0039


Mira3 Quant using 4xnum spec

Spectrum  Sum           O             Si            Pb
K227-5kV  0.9650±0.0136 0.1579±0.0133 0.0925±0.0006 0.7146±0.0025

Actual
          1.0           0.1649        0.0934        0.7426
"""


__revision__ = "$Id: sim-quant-bn-5kV.py John R. Minter $"
__version__ = "0.0.1"


import dtsa2 as dt2
import gov.nist.microanalysis.dtsa2 as gdtsa2
import dtsa2.mcSimulate3 as mc3
import dtsa2.jmGen as jmg
import dtsa2.jmMC3 as jm3

start = time.time()

DataManager.clearSpectrumList()


# change the datDir to match your system
# MacOS
datDir   = "/Users/jrminter/dat/k227-quant"
# Windows10
# datDir =  "D:/dtsa2_sims/dat/k227-quant"

det = findDetector("Oxford p4 05eV 2K")
# det = findDetector("Brucker 5 eV")
e0       =      5.0  # kV
nTraj    =  10000    # trajectories
lt       =    500    # sec
pc       =      5.0  # nA
dose     = pc * lt   # na-sec"
bSaveSpc = True
# coating  = False

# define materials
crocoite = material("PbCrO4", density=6.0)
crocoite.setName("Crocoite")

si = material("Si", density = 2.3290)

k227 = epq.Material(epq.Composition([epq.Element.O,
                                     epq.Element.Si,
                                     epq.Element.Pb],
                                     [ 0.1640,
                                       0.0934,
                                       0.7246 ]),
                                      epq.ToSI.gPerCC(4.0))
k227.setName("k227")

# Start with k227
xrts = []
trs = mc3.suggestTransitions(k227, e0)
for tr in trs:
  xrts.append(tr)

print(trs)

spc_k227 = jm3.simBulkStd(k227, det, e0, nTraj, 100, 1.0, False)
spc_k227.display()
spc_k227.rename("k227-5kV")
spc_k227.setAsStandard(k227)
fi = datDir + "/spc_k227.msa"
print(fi)
if(bSaveSpc):
  spc_k227.save(fi)

# Next Si
xrts = []
trs = mc3.suggestTransitions(si, e0)
for tr in trs:
  xrts.append(tr)

print(trs)

spc_si = jm3.simBulkStd(si, det, e0, nTraj, 100, 1.0, False)
spc_si.display()
spc_si.rename("Si-5kV")
spc_si.setAsStandard(si)
fi = datDir + "/spc_si.msa"
print(fi)
if(bSaveSpc):
  spc_si.save(fi)

# Finally, Crocoite
xrts = []
trs = mc3.suggestTransitions(crocoite, e0)
for tr in trs:
  xrts.append(tr)

print(trs)

spc_crocoite = jm3.simBulkStd(crocoite, det, e0, nTraj, 100, 1.0, False)
spc_crocoite.rename("Crocoite-5kV")
spc_crocoite.setAsStandard(crocoite)
spc_crocoite.display()
fi = datDir + "/spc_crocoite.msa"
print(fi)
if(bSaveSpc):
  spc_crocoite.save(fi)

# Run a quantification - that fails...

qua = quantify(spc_k227, {"Si":spc_si, "Crocoite":spc_crocoite},
               refs={}, preferred=(), elmByDiff=None,
               oByStoic=True, oxidizer=None, extraKRatios=None, fiat={})
tmp = qua.getComposition()
print(tmp.descriptiveString(False))



end = time.time()

delta = end - start

print(delta)



