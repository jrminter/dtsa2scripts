"""
ionizationRange(mat,[e0=20.0],[alg=epq.ElectronRange.KanayaAndOkayama1972])  Ex: ionizationRange(createMaterial(), 30.0)

Tablulates the mean electron range for each element in the sample and each excitable shell.

"""
c = material("C", 2.62)
al = material("Al", 2.70)
fe = material("Fe", 7.86)
cu = material("Cu", 8.96)
ag = material("Ag", 10.5)
au = material("Au", 19.3)
gd = material("Gd", 7.89)

lEl = [c,al,fe,cu,ag,au,gd]
le0 = [3.0, 5.0, 10.0, 20.0, 30.0]

for el in lEl:
	print(el.getName())
	for e0 in le0:
		er = ionizationRange(el, e0)
		print(er)



