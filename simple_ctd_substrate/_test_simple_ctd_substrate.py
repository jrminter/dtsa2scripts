# -*- coding: utf-8 -*-

"""
DTSA-II Script - J. R. Minter

_test_simple_ctd_substrate.py

A DTSA-II Monte Carlo simulation of a simple coated substrate

  Date      Who   Ver    Comment
----------  ---  -----   ------------------------------------------
2020-06-03  JRM  0.0.1   Initial test.

Kelvin:
10-nm-C-on-BN-5-kV = [B(0.5627±0.1133 mass frac),N(0.6475±0.0970 mass frac),Σ=1.2101±0.2103]

NIST DTSA-II Lorentz 2019-07-29:
[OpenJDK 64-Bit Server VM (Amazon.com Inc.)] on java11.0.4

Unmodified
fails with  error: "global name 'coating' is not defined"
Workaround: comment out lines 1266, 1267, 1321
10-nm-C-on-BN-5-kV = [B(0.5689±0.1153 mass frac),N(0.6543±0.0984 mass frac),Σ=1.2232±0.2137]

NIST DTSA-II Lorentz 2020-05-26
Also java11.0.4

Fails with:
Traceback (most recent call last):
  File "<stdin>", line 8, in <module>
  File "/Users/jrminter/Documents/git/dtsa2Scripts/_test_simple_ctd_substrate/_test_simple_ctd_substrate.py", line 95, in <module>
    qua = quantify(spc_bn, {"B":spc_b, "N":spc_aln},
  File "/Applications/NIST DTSA-II Lorentz 2020-05-26/Lib/dtsa2/__init__.py", line 1320, in quantify
    return qus.compute(unknown)
	at Jama.SingularValueDecomposition.<init>(SingularValueDecomposition.java:181)
	at Jama.Matrix.svd(Matrix.java:797)
	at gov.nist.microanalysis.Utility.LinearLeastSquares.perform(LinearLeastSquares.java:87)
	at gov.nist.microanalysis.EPQLibrary.FilterFit.perform(FilterFit.java:270)
	at gov.nist.microanalysis.EPQLibrary.FilterFit.updateUnknown(FilterFit.java:572)
	at gov.nist.microanalysis.EPQLibrary.FilterFit.getKRatios(FilterFit.java:410)
	at gov.nist.microanalysis.EPQLibrary.QuantifyUsingStandards.compute(QuantifyUsingStandards.java:508)
	at java.base/jdk.internal.reflect.NativeMethodAccessorImpl.invoke0(Native Method)
	at java.base/jdk.internal.reflect.NativeMethodAccessorImpl.invoke(NativeMethodAccessorImpl.java:62)
	at java.base/jdk.internal.reflect.DelegatingMethodAccessorImpl.invoke(DelegatingMethodAccessorImpl.java:43)
	at java.base/java.lang.reflect.Method.invoke(Method.java:566)
java.lang.ArrayIndexOutOfBoundsException: java.lang.ArrayIndexOutOfBoundsException: Index -1 out of bounds for length 0
Elapse: 0:04:39.5



Running /Users/jrminter/Documents/git/dtsa2Scripts/_test_simple_ctd_substrate/_test_simple_ctd_substrate.py
10-nm-C-on-BN-5-kV = [B(0.5651±0.1142 mass frac),N(0.6506±0.0976 mass frac),Σ=1.2157±0.2118]
258.210000038
Elapse: 0:04:18.2
"""


__revision__ = "$Id: _test_simple_ctd_substrate.py John R. Minter $"
__version__ = "0.0.1"

import dtsa2 as dt2
import dtsa2.mcSimulate3 as mc
import time

DataManager.clearSpectrumList()

start = time.time()

# MacOS

datDir =  "/Users/jrminter/dat/dtsa2-eds-sim/sim-bn-quant/sim-ctd-substrate/"

det      = findDetector("Bruker 5 eV")
e0       =     2.5  # kV
nTraj    = 10000    # trajectories
lt       =   500    # sec
pc       =     5.0  # nA
dose     = pc * lt  # na-sec"
bSaveSpc = True
tCnm     = 10       # nm

c  = material("C", density=2.1)
b  = material("B", density = 2.37)
bn = material("BN", density = 2.1)
aln = material("AlN", density = 3.26)


"""
Simulate a spectrum of a coating on a substrate.
Monte Carlo simulate a spectrum from a 'substrate' material coated
with 'coating' of the specified thickness (in microns).

coatedSubstrate(coating, thickness, substrate, det, e0=20.0,
                       withPoisson=True, nTraj=defaultNumTraj,
                       dose=defaultDose, sf=defaultCharFluor,
                       bf=defaultBremFluor, xtraParams=defaultXtraParams)

"""

spc_bn = mc.coatedSubstrate(c, 20.0e-03, bn, det, 5.0, True, nTraj, dose, True, True, {})
spc_bn.rename("10-nm-C-on-BN-5-kV")
spc_bn.setAsStandard(bn)
spc_bn.display()

spc_b = mc.coatedSubstrate(c, 20.0e-03, b, det, 5.0, True, nTraj, dose, True, True, {})
spc_b.rename("10-nm-C-on-B-5-kV")
spc_b.setAsStandard(b)
spc_b.display()

spc_aln = mc.coatedSubstrate(c, 20.0e-03, aln, det, 5.0, True, nTraj, dose, True, True, {})
spc_aln.rename("10-nm-C-on-AlN-5-kV")
spc_aln.setAsStandard(aln)
spc_aln.display()

# Run a quantification - that fails...

qua = quantify(spc_bn, {"B":spc_b, "N":spc_aln},
               refs={}, preferred=(), elmByDiff=None,
               oByStoic=False, oxidizer=None, extraKRatios=None, fiat={})
tmp = qua.getComposition()
print(tmp.descriptiveString(False))

end = time.time()
delta = end - start

print(delta)






