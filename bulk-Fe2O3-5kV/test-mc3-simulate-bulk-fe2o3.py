# test-mc3-simulate-fe2o3
#
# 2021-06-10 J.R. Minter
#
import os, shutil
import dtsa2.mcSimulate3 as mc3

start = time.time()

homDir   = os.environ['HOME']
homDir   = homDir.replace('\\','/')
# store the results in the outDir directory (a place for spectra...)
outDir   = homDir + "/Documents/git/dtsa2Scripts/utility/output/"

bSaveSpc = True
DataManager.clearSpectrumList()


e0    = 5.0
nTraj = 10000
det   = findDetector("Oxford p4 05eV 4K")
mat   = material("Fe2O3",5.0)

spc = mc3.simulate(mat, det, e0, 120.0, True, nTraj, True, True, {})
sName = "Bulk Fe2O3"
spc.rename(sName)
spc.display()
fi = outDir + "/Fe2O3-%g-kV-%g-Traj.msa" % (e0, nTraj)
spc.save(fi)


end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg