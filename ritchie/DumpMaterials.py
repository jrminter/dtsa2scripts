# -*- coding: utf-8 -*-
# A script for dumping all the standards currently defined in the database to a Python script.
# The Python script can be used to export the standards to DTSA-II on another computer.
# Nicholas W. M. Ritchie 21-Feb-2011

# Change this to point to where you want the resulting script file to be saved.
filename = jl.System.getProperty("user.home")+"/stds.py"

# Get a set containing all standards
stds=Database.findAllStandards()

# Build up the script in the variable 'text'
text = "# -*- coding: utf-8 -*-\n\n"
# Define a helper method 'defineStd'
text += "def defineStd(elms,qty,name,density=None):\n"
text += "\tc=epq.Composition(map(element,elms),qty,name)\n"
text += "\tif density:\n"
text += "\t\tc=epq.Material(c,epq.ToSI.gPerCC(density))\n"
text += "\tDatabase.addStandard(c)\n\n"
# For each standard in the database
for mat in stds:
	elms = "(%s,)" % (",".join(["\"%s\"" % elm.toAbbrev() for elm in mat.getElementSet()]))
	qty = "(%s,)" % (",".join(["%f" % mat.weightFraction(elm,False) for elm in mat.getElementSet()]))
	text += "defineStd(%s,%s,\"" % (elms, qty)
	text += mat.toString()
	text += "\""
	if isinstance(mat,epq.Material):
	   text += ",%f)\n" % epq.FromSI.gPerCC(mat.getDensity())
	else:
		text += ")\n"

# Write the result to a text file.
print text
import codecs
f = codecs.open(filename, 'w', encoding='utf-8')
f.write(text)
f.close()

