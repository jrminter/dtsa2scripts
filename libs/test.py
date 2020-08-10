import dtsa2 as dtsa2
import dtsa2.jmGen as jmg
e0 = 5.0
si = material("Si", density=2.329)
sio2 = material("SiO2", 2.65)
det = findDetector("Si (Li)")
spc = simCoatedSubstrate(si, sio2, 20, det, e0, 100, "C:/Users/johnr/Documents", 5.0e-6, 100, 1.0, 512)
spc.display()