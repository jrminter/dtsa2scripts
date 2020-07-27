# -*- coding: utf-8 -*-
__revision__ = "$Id: jrmFoo.py JRM $"
__version__ = "0.0.2 - 2017-05-26"

"""

2015-05-25 on crunch
2015-05-25 on win10vm too

Jython 2.7.0 (default:9987c746f838, Apr 29 2015, 02:25:11) 
[Java HotSpot(TM) 64-Bit Server VM (Oracle Corporation)] on java1.8.0_131
Welcome to DTSA-II's scripting interface.  Type 'help()' for more information.
<Ctrl-Return> executes a command.  <Return> to create multi-line commands.
Elapse: 0:00:00.1
1> import dtsa2.jrmFoo as jf
2> dir(jf)
['__builtins__', '__doc__', '__file__', '__name__', '__package__',
 '__revision__', '__version__', 'dtsa2', 'epd', 'epq', 'ept', 'epu',
 'iio', 'jio', 'jl', 'nm', 'nm3', 'rp', 'sys']
3> jf.rp()
u'C:\\Users\\jrminter\\Documents\\DTSA\\2017\\May\\25-May-2017'
4> 

"""

import sys
sys.packageManager.makeJavaPackage("gov.nist.microanalysis.NISTMonte.Gen3", "CharacteristicXRayGeneration3, BremsstrahlungXRayGeneration3, FluorescenceXRayGeneration3, XRayTransport3", None)
import gov.nist.microanalysis.NISTMonte as nm
import gov.nist.microanalysis.NISTMonte.Gen3 as nm3
import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQLibrary.Detector as epd
import gov.nist.microanalysis.Utility as epu
import gov.nist.microanalysis.EPQTools as ept
import javax.imageio.ImageIO as iio
import java.io as jio
import java.lang as jl
import dtsa2

def rp():
   """run(script)
   Runs the script in the file specified by script."""
   a = dtsa2.reportPath()
   return a


