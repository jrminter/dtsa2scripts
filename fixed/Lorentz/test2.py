import dtsa2 as dtsa2
import dtsa2.jmMC3 as jm3
import dtsa2.jmGen as jmg
import time

DataManager.clearSpectrumList()



start = jmg.startTimer()

outPath = "C:/Users/johnr/Documents/test"
cu    = material("Cu", density=8.92)
e0    = 15.0
nTraj = 10000
det = findDetector("Oxford p4 05eV 4K")
#                                                       im dim  lt  pc emiSz
cu_spec = jm3.uncoatedSimBulkStd(cu,det,e0,nTraj,outPath,5.0e-6,100,1.0,512)
cu_spec.display()
cu_spec.rename("Cu 15kV")
fi =  outPath + "/"
fi += "%s-%g-kV.msa" % ("Cu", e0)
cu_spec.save(fi)

jmg.endTimer(start)
