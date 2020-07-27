# -*- coding: utf-8 -*-

# Simulate standards and thin film specimens of ZnO on Si for
# measuring K-ratios for thin film measurement.
#
# For 500 trajectories at 4 kV:
# This script required 1.248 min
# Elapse: 0:01:14.9
#
# For 50000 trajectories at 4kV:
# This script required 144.035 min
# ...or 2.401 hr
# Elapse: 2:24:02.7
import sys
import os
import time
import shutil
import gov.nist.microanalysis.Utility as nu
import dtsa2.mcSimulate3 as mc3
import dtsa2.jmGen as jmg

# adjustable parameters
e0      = 4.0
det     = findDetector("Oxford p4 05eV 2K")
lThNm   = [  0.1,   5.0,  10.0,  20.0,  30.0,  40.0,  50.0, 60.0,
            70.0,  80.0,  90.0, 100.0, 110.0, 115.0, 120.0,
           125.0, 130.0, 135.0, 140.0, 145.0, 150.0]
nTraj   = 50000     # trajectories
dose    =  5000   # nA-sec
nDigits = 6


# directories
gitDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/utility"
prjDir = gitDir + relPrj
resDir = prjDir + '/sim-zno-on-si-4kV'
jmg.ensureDir(resDir)
rptDir = prjDir + '/sim-zno-on-si-4kV Results/'


DataManager.clearSpectrumList()
start = time.time()

zno = material("ZnO", density=5.61)
si  = material("Si", density=2.33)

xrts = []
trs = mc3.suggestTransitions(zno, e0)
for tr in trs:
    xrts.append(tr)

xtraParams={}
xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))
# xtraParams.update(mc3.configureOutput(simDir))
spc_zno_std = mc3.simulate(zno, det, e0, dose, True, nTraj, True, True, xtraParams)
sName = "ZnO std-%g-kV" % (e0)
spc_zno_std.rename(sName)
spc_zno_std.setAsStandard(zno)
spc_zno_std.display()
fi =  resDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (nTraj)
spc_zno_std.save(fi)

xrts = []
trs = mc3.suggestTransitions(si, e0)
for tr in trs:
    xrts.append(tr)

xtraParams={}
xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))
# xtraParams.update(mc3.configureOutput(simDir))
spc_si_std = mc3.simulate(si, det, e0, dose, True, nTraj, True, True, xtraParams)
sName = "Si std-%g-kV" % (e0)
spc_si_std.rename(sName)
spc_si_std.setAsStandard(si)
spc_si_std.display()
fi =  resDir + "/"
fi += sName
fi += "-%g-Traj.msa" % (nTraj)
spc_si_std.save(fi)

xrts = []
trs = mc3.suggestTransitions(zno, e0)
for tr in trs:
    xrts.append(tr)

trs = mc3.suggestTransitions(si, e0)
for tr in trs:
    xrts.append(tr)

xtraParams={}
xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))

stds =  { element("O"): spc_zno_std, element("Zn"): spc_zno_std, element("Si"): spc_si_std }

print("Specimen Name \t\t Oxygen Ka1 K-Ratio \t Zinc La1 K-Ratio \t Si Ka K-Ratio ")
      #5 nm ZnO on Si at 4 kV  0.073166 ± 0.000133 0.075862 ± 0.000106 0.877289 ± 0.000202

for thNm in lThNm:
    thUm = thNm/1000.
    spc = mc3.coatedSubstrate(zno, thUm, si, det, e0,
                              True, nTraj, dose, True,
                              True, xtraParams)
    sName = "%g nm ZnO on Si at %g kV" % (thNm, e0)
    spc.rename(sName)
    spc.display()
    fi =  resDir + "/"
    fi += sName
    fi += "-%g-Traj.msa" % (nTraj)
    spc.save(fi)

    ff=epq.FilterFit(det,epq.ToSI.keV(e0))
    ff.addReference(element("O"),  spc_zno_std)
    ff.addReference(element("Zn"), spc_zno_std)
    ff.addReference(element("Si"), spc_si_std)
    kr  = ff.getKRatios(spc)
    kO  = kr.getKRatioU(epq.XRayTransition(epq.Element.O, epq.XRayTransition.KA1))
    kZn = kr.getKRatioU(epq.XRayTransition(epq.Element.Zn,epq.XRayTransition.LA1))
    kSi = kr.getKRatioU(epq.XRayTransition(epq.Element.Si,epq.XRayTransition.KA1))
    print (u"%s\t%g \u00B1 %g\t%g \u00B1 %g\t%g \u00B1 %g" % ( sName,
                                                               round(kO.doubleValue(), nDigits ),
                                                               round(kO.uncertainty(), nDigits ),
                                                               round(kZn.doubleValue(), nDigits ),
                                                               round(kZn.uncertainty(), nDigits ),
                                                               round(kSi.doubleValue(), nDigits ),
                                                               round(kSi.uncertainty(), nDigits )))

    # kr = kratios(spc, stds, refs={})
    # print(kr[0])




# clean up cruft
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



