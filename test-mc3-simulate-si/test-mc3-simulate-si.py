# C:\Users\johnr\Documents\git\dtsa2scripts\test-mc3-simulate-si\test-mc3-simulate-si.py
#
# 2021-06-11 J.R. Minter (updated)
#
# at 5kV with 10000 trajectories: 
#
import os, shutil, time
import dtsa2.mcSimulate3 as mc3

start = time.time()

homDir = os.environ['HOME']
homDir = homDir.replace('\\','/')
wrkDir = homDir + "/Documents/git/dtsa2Scripts/utility"
outDir = homDir + "/Documents/git/dtsa2Scripts/utility/output/"
rptDir = wrkDir + '/test-mc3-simulate-si/'

# Note: On Win10 MSA files are in:
# C:\Users\johnr\Documents\git\dtsa2scripts\utility\output
# 10000 trajectories at  5 kV This script required 1.244 min
# 10000 trajectories at 20 kV This script required 4.205 min


bSaveSpc = True
DataManager.clearSpectrumList()

e0    = 20.0  # or 5.0
nTraj = 10000
det   = findDetector("Oxford p4 05eV 4K")
mat   = material("Si", 2.329)

spc = mc3.simulate(mat, det, e0, 120.0, True, nTraj, True, True, {})
sName = "Bulk-Si-%g" % (e0)
spc.rename(sName)
spc.display()
fi = outDir + "/Si-%g-kV-%g-Traj.msa" % (e0, nTraj)
spc.save(fi)

print "Done!"

end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg