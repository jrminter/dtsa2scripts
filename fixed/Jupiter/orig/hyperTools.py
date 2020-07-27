"""Tools for processing hyperspectral data sets."""

import gov.nist.microanalysis.EPQLibrary as epq
import gov.nist.microanalysis.EPQTools as ept
import dtsa2 as dt2
import java.lang as jl

def openRipple(rpl, e0, i0, liveTime, det):
   """openRipple(rpl, e0, i0, liveTime, det)
   Open a ripple file as an ISpectrumData object with setPosition(x,y) to permit navigating /
through the spectra."""
   sp = epq.SpectrumProperties()
   sp.setDetector(det)
   sp.setNumericProperty(epq.SpectrumProperties.BeamEnergy, e0)
   sp.setNumericProperty(epq.SpectrumProperties.FaradayBegin, i0)
   sp.setNumericProperty(epq.SpectrumProperties.LiveTime, liveTime)
   return ept.RippleSpectrum(rpl, sp)

def maskRipple(inRpl, outFile, mask):
   """maskRipple(inRpl, outFile, mask)
   Sets the individual data items to zero based on the specified mask.  If mask.getRGB(c,r)>0 /
then copy the contents at(c,r) of inRpl to outFile.rpl.  Otherwise the contents of outFile /
is set to all zeros."""
   outRpl = "%s.rpl" % outFile
   outRaw = "%s.raw" % outFile
   len = rpl.getDepth()
   ty = rpl.getDataType()
   res = ept.RippleFile(rpl.getColumns(), rpl.getRows(), rpl.getDepth(), rpl.getDataType(), rpl.getDataSize(), ept.RippleFile.DONT_CARE_ENDIAN, outRpl, outRaw)
   zero = (0) * len
   for c in xrange(0, rpl.getColumns()):
      for r in xrange(0, rpl.getRows()):
         rpl.setPosition(c, r)
         res.setPosition(c, r)
         if mask.getRGB(c, r) > 0:
            if ty == rpl.FLOAT:
               res.write(rpl.readDouble(len))
            else:
               res.write(rpl.readInt(len))
   return res

def quantify(rpl, stds, refs={}, preferred=(), elmByDiff=None, elmByStoic=None, assumedStoic={}, mask=None, step=1, visualize=False, zaf=None, withUnc=False):
   """quantify(rpl,stds,[refs={}],[preferred=()],[elmByDiff=None],[elmByStoic=None],[assumedStoic={}], [mask=None],[zaf=None], [withUnc=False])
   Quantify a ripple/raw spectrum object based on the standards, references and other parameters specified. /
The arguments are the same as dtsa2.multiQuant.  An additional 'mask' argument allows uninteresting pixels /
to be ignored (not quantified.)  The mask should be an object like a BufferedImage with a getRGB(x,y) method. /
The pixel is ignored if mask.getRGB(x,y)==0. The result is written to a RPL/RAW file in FLOAT format.
> import javax.imageio.ImageIO as io
> mask = io.read(jio.File("c:/image.png"))
> zaf = epq.CorrectionAlgorithm.NullCorrection for k-ratios"""
   oldSt = None
   try:
      if (zaf != None) and isinstance(zaf, epq.CorrectionAlgorithm):
         oldSt = epq.AlgorithmUser.getGlobalStrategy()
         newSt = epq.AlgorithmUser.getGlobalStrategy()
         newSt.addAlgorithm(epq.CorrectionAlgorithm, zaf)
         epq.AlgorithmUser.applyGlobalOverride(newSt)         
      det = rpl.getProperties().getDetector()
      e0 = rpl.getProperties().getNumericProperty(epq.SpectrumProperties.BeamEnergy)
      mq = dt2.multiQuant(det, e0, stds, refs, preferred, elmByDiff, elmByStoic, assumedStoic)
      base = rpl.getProperties().getTextProperty(epq.SpectrumProperties.SourceFile)
      compRpl = base.replace(".rpl", "_comp.rpl")
      compRaw = base.replace(".rpl", "_comp.raw")
      compTxt = base.replace(".rpl", "_comp.txt")
      status = file(compTxt, "wt")
      status.write("File:\t%s\n" % base)
      status.write("Results:\t%s\n" % compRpl)
      status.write("Detector:\t%s\n" % det)
      status.write("Beam energy\t%g keV\n" % e0)
      status.write("Standards\n")
      i = 0;
      elms = []  # Ensures the plane order is correct
      if det.getChannelCount() != rpl.getChannelCount():
            print "ERROR: The number of channels in %s (%d) doesn't match the number of channels in the RPL file (%d)." % (det, det.getChannelCount(), rpl.getChannelCount())
            return 
      for elm, std in stds.iteritems():
         if std.getChannelCount() != rpl.getChannelCount():
            print "ERROR: The number of channels in %s (%d) doesn't match the number of channels in the RPL file (%d)." % (std, std.getChannelCount(), rpl.getChannelCount())
            return 
         status.write("\t%d\t%s\t%s\n" % (i, elm, std))
         elms.append(dt2.element(elm))
         i = i + 1
      if len(refs) > 0:
         status.write("References\n")
         for xrt, ref in refs.iteritems():
            if ref.getChannelCount() != rpl.getChannelCount():
               print "ERROR: The number of channels in %s (%d) doesn't match the number of channels in the RPL file (%d)." % (ref, ref.getChannelCount(), rpl.getChannelCount())
            status.write("\t%s\t%s\n" % (xrt, ref))
      if len(preferred) > 0:
         status.write("Preferred transitions\n")
         for xrt in preferred:
            status.write("\t%s for %s\n" % (xrt, xrt.getElement()))
      if elmByDiff:
         status.write("Element by difference: %s\n" % elmByDiff)
      if elmByStoic:
         status.write("Element by Stoiciometry: %s\n" % elmByStoic)
         status.write("Element\tValence\n")
         for elm, stoic in assumedStoic.iteritems():
            status.write("\t%s\t%g\n" % (elm, stoic))
      comps = ept.RippleFile((rpl.getColumns() + step - 1) / step, (rpl.getRows() + step - 1) / step, len(stds) + 1, ept.RippleFile.FLOAT, 8, ept.RippleFile.DONT_CARE_ENDIAN, compRpl, compRaw)
      uncert = None
      if withUnc:
          uncRaw = base.replace(".rpl", "_unc.raw")
          uncRpl = base.replace(".rpl", "_unc.rpl")
          uncert = ept.RippleFile((rpl.getColumns() + step - 1) / step, (rpl.getRows() + step - 1) / step, len(stds) + 1, ept.RippleFile.FLOAT, 8, ept.RippleFile.DONT_CARE_ENDIAN, uncRpl, uncRaw)
      dumpIt = False
      if dumpIt:
         dumpRpl = base.replace(".rpl", "_dump.rpl")
         dumpRaw = base.replace(".rpl", "_dump.raw")
         dump = ept.RippleFile((rpl.getColumns() + step - 1) / step, (rpl.getRows() + step - 1) / step, rpl.getChannelCount(), ept.RippleFile.UNSIGNED, 4, ept.RippleFile.DONT_CARE_ENDIAN, dumpRpl, dumpRaw)
      rpl.setSpan(step, step)
      for r in xrange(0, rpl.getRows(), step):
         if dt2.isTerminated():
            break
         dt2.StdOut.append(".")
         dt2.StdOut.flush()
         for c in xrange(0, rpl.getColumns(), step):
            if dt2.isTerminated():
               break
            if visualize:
               dt2.clearSpectra()
               dt2.display(rpl)
            comps.setPosition(c / step, r / step)
            if (mask == None) or ((mask.getRGB(c, r) & 0xFFFFFF) > 0):
                try:
                   rpl.setPosition(c, r)
                   rs = epq.SpectrumUtils.copy(rpl)
                   if dumpIt:
                      # print "%d\t%d\t%d\t%d" % (c, r, c / step, r / step)
                      dump.setPosition(c / step, r / step)
                      dump.write(epq.SpectrumUtils.toIntArray(rs))
                   res = mq.compute(rs)
                   comp = res.getComposition()
                   if visualize:
                      rpl.getProperties().setCompositionProperty(epq.SpectrumProperties.MicroanalyticalComposition, comp)
                      dt2.annotComposition()
                   tmp, unc = [], []
                   sU = 0.0
                   for elm in elms:
                      tmp.append(comp.weightFraction(elm, True))
                      u = comp.weightFractionU(elm, True).uncertainty()
                      unc.append(u)
                      sU = sU + u*u
                   tmp.append(comp.sumWeightFraction())
                   unc.append(jl.Math.sqrt(sU))
                   if dumpIt:
                      print tmp
                   comps.write(tmp)
                   if uncert:
                       uncert.write(unc)
                except (epq.EPQException, Exception, jl.Exception), e:
                   msg = "row = %d, col = %d failed: %s" % (r, c, e)
                   print msg
                   status.write(msg + "\n")
                   for elm, std in stds.iteritems():
                      comps.write(0.0)
                      if uncert:
                         uncert.write(unc)
                   if visualize:
                      dt2.setAnnotation("row = %d, col = %d failed." % (r, c))
            else:
               for elm, std in stds.iteritems():
                  comps.write(0.0)
                  if uncert:
                     uncert.write(unc)
               comps.write(1.0)
               if uncert:
                  uncert.write(1.0)
      comps.close()
      if uncert:
         uncert.close()
      if dumpIt:
         dump.close()
      status.close()
   finally:
      if oldSt:
         epq.AlgorithmUser.applyGlobalOverride(oldSt)
   
   print "\nDone!"
   
def maxPixel(rpl):
   """maxPixel(rpl)
   Computes the max pixel spectrum for the specified ripple/raw spectrum object."""
   xs = epq.ExtremumSpectrum()
   for r in xrange(0, rpl.getRows()):
      dt2.StdOut.append(".")
      if dt2.terminated:
         break
      for c in xrange(0, rpl.getColumns()):
         rpl.setPosition(r, c)
         xs.include(rpl)
   return xs

def rippleToTIFF(rpl, minLayer=0, maxLayer=4294967295):
   """rippleToTIFF(rpl, min=0, max=4294967295)
   Converts a RPL/RAW file pair to a TIFF image representing the sum of layers minLayer to maxLayer.  
   The rpl argument is the full path of the RPL file.  min and max bound the layers to export."""
   rf = ept.RippleFile(rpl, True)
   tiff = rpl.replace(".rpl", "-[%d,%d].tif" % (minLayer, maxLayer,))
   rf.layersToTIFF(tiff, minLayer, maxLayer+1, True)

      
def rippleToMaps(rpl, vecSet, step=1):
    """RippleToMaps(rpl, vecSet, step=1)
    Applies a vector set to a RPL/RAW file pair producing a dictionary of MapImage objects /
indexed by RegionOfInterest objects."""
    width = rpl.getColumns()
    height = rpl.getRows()
    vecs = vecSet.getVectors()
    planes = epq.MapImage(width, height, len(vecs), str(vecSet))
    rpl.setSpan(step, step)
    for r in xrange(0, height, step):
        if dt2.isTerminated():
            break
        for c in xrange(0, width, step):
            if dt2.isTerminated():
                break
            rpl.setPosition(r, c)
            krs = vecSet.getKRatios(rpl)
            kra = []
            for vec in vecs:
                kra.append(krs.getKRatio(vec.getROI()))
            planes.inc(c,r,kra)
    return planes
