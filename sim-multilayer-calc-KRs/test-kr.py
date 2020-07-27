# -*- coding: utf-8 -*-
#         1         2         3         4         5         6         7 |
# 23456789012345678901234567890123456789012345678901234567890123456789012
#
# test-kr.py
# 2018-09-24 - move KR calc to 
#              jmg.compKRs(spc, stds, trs, det, e0, digits=5)

import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3, FluorescenceXRayGeneration3, XRayTransport3", None)
import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQLibrary.Detector as epd
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.microanalysis.NISTMonte.Gen3 as nm3
import gov.nist.microanalysis.EPQTools as et
import dtsa2.mcSimulate3 as mc3
import dtsa2.jmGen as jmg
import dtsa2.jmMC3 as jm3
import java.util as jutil
import java.io as jio
import java.nio.charset as cs
import os
import shutil
import time
import datetime

det = findDetector("Oxford p4 05eV 2K")
e0       =     4    # kV
nTraj    = 50000    # trajectories
lt       =   500    # sec
pc       =     5.0 # nA
dose     = pc * lt  # na-sec"
bSaveSpc = True
digits   = 5

# To process all data
#
lNmCsim = [  1.0,   5.0,  10.0,  15.0,  20.0,  25.0,  30.0,  35.0,  40.0,
            45.0,  50.0,  60.0,  70.0,  80.0,  90.0, 100.0, 110.0, 120.0, 
           130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0]

# For testing
# start with 20 nm
# lNmCsim = [20.0]

start = time.time()

DataManager.clearSpectrumList()

homDir = os.environ['HOME']
homDir = homDir.replace('\\', '/')
relPrj = "/Documents/work/dtsa2/sim-C-on-Fe-4kV"
datDir = homDir + relPrj + "/msa-%g" % (nTraj)
jmg.ensureDir(datDir)
rptDir = homDir + relPrj + "/test-kr Results"
csvFil = homDir + relPrj + "/dtsa2-C-on-Fe-%g-kV-kratios-%g-traj.csv" % (e0, nTraj)

sName = "%s-std-%g-kV-%g-traj" % ("Fe", e0, nTraj)
sPath = datDir + "/" + sName +".msa"
print(sPath)
# sFile = os.path.join(path, file)
fe_std_spc = readSpectrum(sPath)
display(fe_std_spc)

sName = "%s-std-%g-kV-%g-traj" % ("C", e0, nTraj)
sPath = datDir + "/" + sName +".msa"
# print(sPath)
# sFile = os.path.join(path, file)
c_std_spc = readSpectrum(sPath)
display(c_std_spc)


sName = "%s-std-%g-kV-%g-traj" % ("Fe", e0, nTraj)
sPath = datDir + "/" + sName +".msa"
# print(sPath)
# sFile = os.path.join(path, file)
fe_std_spc = readSpectrum(sPath)
display(fe_std_spc)

# only get FeL at 7 kV
trs = [epq.XRayTransitionSet(epq.Element.C, epq.XRayTransitionSet.K_FAMILY),
       epq.XRayTransitionSet(epq.Element.Fe, epq.XRayTransitionSet.L_FAMILY)]
cStd  = {"El":element("C"),  "Spc":c_std_spc}
feStd = {"El":element("Fe"), "Spc":fe_std_spc}

stds = [cStd, feStd]

nSpec = len(lNmCsim)

l_nm_C = []
l_kC_mu  = []
l_kC_unc  = []
l_kFe_mu = []
l_kFe_unc = []

for tc in lNmCsim:
    sName = "%g-nm-C-on-Fe-%g-kV-%g-traj" % (tc, e0, nTraj)
    fi =  datDir + "/"
    fi += sName
    fi += ".msa"
    # print(fi)
    spc = readSpectrum(fi)
    display(spc)
    kr = jmg.compKRs(spc, stds, trs, det, e0, digits=5)
    l_nm_C.append(tc)
    l_kC_mu.append(kr[0][0])
    l_kC_unc.append(kr[0][1])
    l_kFe_mu.append(kr[1][0])
    l_kFe_unc.append(kr[1][1])
    # print(kr)

f=open(csvFil, 'w')
strLine = 'nm_c,c_ka_mu,c_ka_unc,fe_la_mu,fe_la_unc\n'
f.write(strLine)
for i in range(nSpec):
    strLine = "%g" % l_nm_C[i] + ","
    strLine = strLine + "%.5f" % l_kC_mu[i] + ","
    strLine = strLine + "%.5f" % l_kC_unc[i] + ","
    strLine = strLine + "%.5f" % l_kFe_mu[i] + ","
    strLine = strLine + "%.5f" % l_kFe_unc[i] + "\n"
    f.write(strLine)  
f.close()

# clean up cruft
shutil.rmtree(rptDir)
print "Done!"

end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg

