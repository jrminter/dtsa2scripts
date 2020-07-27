import dtsa2 as dtsa2
import dtsa2.mcSimulate3 as mc3
import dtsa2.jmMC3 

DataManager.clearSpectrumList()


outPath = "C:/Users/johnr/Documents/git/dtsa2Scripts/ben-buse/out"

test = {}

# GMR Film does not like 20 kV - right on a limit for Ag, so
# let's go to 15 kV
e0  = 15 # accelerating voltage

# det = findDetector("Si(Li)") # Replace with your detector's name
det   = findDetector("Oxford p4 05eV 4K") # Replace with your detector's name
nTraj = 10000	# Number of electrons simulated

rhoAg = 10.49

# create materials
ag    = epq.Material(epq.Composition(map(dtsa2.element,["Ag"],),[1],"Silver"),
        epq.ToSI.gPerCC(10.5))   # ToSi converts density g/cc to kg/m3, \
                                 # composition as mass fraction
s316H = epq.Material(epq.Composition(map(dtsa2.element,["Mn","Cr","Mo","Ni","Fe"],),
	                                                   [0.016,0.169,0.022,0.119,0.667],
	                                                   "316H steel"),epq.ToSI.gPerCC(7.8))

# calculate electron range for one of the materials - determining range of x-ray emission images and phi-rho-z
range = dtsa2.electronRange(s316H, e0, density=epq.ToSI.gPerCC(7.8))

print(range)

# Start with metal standards
cr = material("Cr", density=7.19)
mn = material("Mn", density=7.43)
ni = material("Ni", density=8.90)
fe = material("Fe", density=7.86)


def fullSimBulkStd(mat, det, e0, nTraj, outPath, dim=5.0e-6, lt=100,
	               pc=1.0, emiSize=512, ctd=False):
    """
    fullSimBulkStd(mat, det, e0, nTraj, outPath, dim=5.0e-6, lt=100,
                   pc=1.0, emiSize=512, ctd=False)

    Use mc3 simulation to simulate an uncoated standard specimen

    Parameters
    ----------
    mat - a dtsa material.
        Note the material must have an associated density. It should
        have a useful name.

    det - a dtsa detector
        Here is where we get the detector properties and calibration

    e0 - float
        The accelerating voltage in kV

    nTraj - integer
        The number of trajectories to run

    outPath - string
        The path to the directory for output

    dim - float (5.0e-6)
        The size of the emission images

    lt - integer (100)
        The live time (sec)

    pc - float (1.0)
        The probe current in nA

    emiSize - int (default 512)
        The width and depth of the emission images.

    ctd - Boolean (False) - is C coated


    Returns
    -------
    sim - DTSA scriptable spectrum 
        The simulated standard spectrum
    
    Example
    -------
    import dtsa2 as dtsa2
    import dtsa2.mcSimulate3 as mc3
    det = findDetector("Oxford p4 05eV 4K")
    cu = material("Cu", density=8.92)

    outPath = "C:/Users/johnr/Documents/git/dtsa2Scripts/ben-buse/out"

    a = fullSimBulkStd(cu, det, 15.0, 100, outPath,
                       dim=5.0e-6, lt=100,
                       pc=1.0, emiSize=512, ctd=False)
    a.display()

    """
    dose = pc * lt  # na-sec"
    xrts = []

    trs = mc3.suggestTransitions(mat, e0)
    for tr in trs:
        xrts.append(tr)
    mc3.configureEmissionImages(xrts, dim, emiSize)

    xtraParams={}
    xtraParams.update(mc3.configureXRayAccumulators(xrts,True, True, True))
    # note that the image size on the specimen is in meters...
    xtraParams.update(mc3.configureEmissionImages(xrts, 5.0e-6, 512))
    xtraParams.update(mc3.configurePhiRhoZ(5.0e-6))
    xtraParams.update(mc3.configureTrajectoryImage(5.0e-6, 512))
    xtraParams.update(mc3.configureVRML(nElectrons=100))
    xtraParams.update(mc3.configureOutput(outPath))
    mc3.configureOutput(outPath)
    print("Output sent to %s") % (outPath)

    dose = lt*pc

    sim = mc3.simulate(mat, det, e0, dose, withPoisson=True, nTraj=nTraj,
                       sf=True, bf=True, xtraParams=xtraParams)

    sName = "%s-%g-kV" % (mat, e0)
    sim.rename(sName)
    sim.setAsStandard(mat)
    sim.display()
    fi =  outPath + "/"
    fi += sName
    fi += "-%g-Traj.msa" % (nTraj)
    print(fi)
    sim.save(fi)
    return(sim)


#             simBulkStd(mat, det, e0, nTraj, lt=100, pc=1.0, ctd=True)
# get 500 combines

spc_cr = 




spc_cr = dtsa2.jmMC3.simBulkStd(cr, det, e0, nE, 100, 5.0)
sName = "Cr-%g" % (e0)
spc_cr.rename(sName)
spc_cr.setAsStandard(cr)
spc_cr.display()
fi =  outPath + "/"
fi += sName
fi += "-%g-Traj.msa" % (nE)
print(fi)
spc_cr.save(fi)


spc_mn = dtsa2.jmMC3.simBulkStd(mn, det, e0, nE, 100, 5.0)
sName = "Mn-%g" % (e0)
spc_mn.rename(sName)
spc_mn.setAsStandard(mn)
spc_mn.display()
fi =  outPath + "/"
fi += sName
fi += "-%g-Traj.msa" % (nE)
spc_mn.save(fi)

spc_fe = dtsa2.jmMC3.simBulkStd(fe, det, e0, nE, 100, 5.0)
sName = "Fe-%g" % (e0)
spc_fe.rename(sName)
spc_fe.setAsStandard(fe)
spc_fe.display()
fi =  outPath + "/"
fi += sName
fi += "-%g-Traj.msa" % (nE)
spc_fe.save(fi)

spc_ni = dtsa2.jmMC3.simBulkStd(ni, det, e0, nE, 100, 5.0)
sName = "Ni-%g" % (e0)
spc_ni.rename(sName)
spc_ni.setAsStandard(fe)
spc_ni.display()
fi =  outPath + "/"
fi += sName
fi += "-%g-Traj.msa" % (nE)
spc_fe.save(fi)


ag = material("Ag", density=10.49)
#             simBulkStd(mat, det, e0, nTraj, lt=100, pc=1.0, ctd=True)
# get 500 combines
spc_ag = dtsa2.jmMC3.simBulkStd(ag, det, e0, nE, 100, 5.0)
sName = "Ag-%g" % (e0)
spc_ag.rename(sName)
spc_ag.setAsStandard(ag)
spc_ag.display()
fi =  outPath + "/"
fi += sName
fi += "-%g-Traj.msa" % (nE)
spc_ag.save(fi)

spc_ss = dtsa2.jmMC3.simBulkStd(s316H, det, e0, nE, 100, 5.0)
sName = "Stainless-Steel-%g" % (e0)
spc_ss.rename(sName)
spc_ss.setAsStandard(s316H)
spc_ss.display()
fi =  outPath + "/"
fi += sName
fi += "-%g-Traj.msa" % (nE)
spc_ss.save(fi)


# set xray transitions
trs = [epq.XRayTransition(epq.Element.Ag, epq.XRayTransition.LB1),
epq.XRayTransition(epq.Element.Fe, epq.XRayTransition.KB1)]

print(trs)


# create samples consisting of silver film on steel (316H) substrate
film = {}
film[1] = [ag , 0.000000020],[s316H, 0.000010]		# 0.000005000 = 5 um 0.000010 = 10um	Configuring multi-layers composition and thickness
film[2] = [ag , 0.000000015],[s316H, 0.000010]		# 0.000000050 = 0.05 um
film[3] = [ag , 0.000000010],[s316H, 0.000010]

xtraP = {}
xtraP = {"Characteristic Accumulator":True, "Char Fluor Accumulator":True, "Brem Fluor Accumulator":True}
# outpathoutPath  = "/Users/jrminter/Documents/git/dtsa2Scripts/ben-buse/out/"
xtraP.update(mc3.configureOutput(outPath))
xtraP.update(mc3.configurePhiRhoZ(1.5*range))
xtraP.update(mc3.configureEmissionImages(trs, 1.5*range, size = 512))
xtraP.update(mc3.configureTrajectoryImage(1.5*range, size = 512))




resF = {}
for fil in film:
	extension = str(film[fil][0])
	extension = extension.replace("[","")
	extension = extension.replace("]","")
	extension = extension.replace(",","")
	# pathloc = 'O:\Documents\PFE_Data\Users\Charles_Younes\\041115_CarbonInSteel\\dtsa2repeat7_' + extension # Change to output folder

	xtraP.update(mc3.configureOutput(outPath)) 	
	print xtraP
	resF[fil] = mc3.multiFilm(film[fil], det,e0=e0, nTraj=nE, dose=500.0, sf=True, bf=True,xtraParams=xtraP) # run simulations
	tmp = str(resF[fil])
	tmpx = extension + tmp
	resF[fil].rename(tmpx)
	resF[fil].save("%s/%s.msa" % ( outPath, resF[fil] ))	# change outPath above to output folder ( location, name) it adds extension
	resF[fil].display()

"""
