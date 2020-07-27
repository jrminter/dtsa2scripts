# -*- coding: utf-8 -*-

# test uncertain value

import sys
import os
import time
import shutil
import gov.nist.microanalysis.Utility as epu
import dtsa2.jmGen as jmg

gitDir = os.environ['GIT_HOME']
relPrj = "/dtsa2Scripts/utility"
prjDir = gitDir + relPrj
rptDir = prjDir + '/test-uv Results/'


start = time.time()

nmZnO1   = 40.1
uvOKa1   = epu.UncertainValue2(0.269157,0.000126)
uvZnLa1  = epu.UncertainValue2(0.259251,9.4e-05)
uvSiKa1  = epu.UncertainValue2(0.654561,8.4e-05)

def make_str_lst_unc_val(id, luv):
    """
    make_str_lst_unc_val(id, luv)

    Make a formatted string from an ID string and a list of uncertain values.

    Input
    -----
    id  A number or a string that will be output as a string.

    luv A list of DTSA-II UncertainValue2 items. These will be printed
        as comma-delimited pairs with 6 digits following the decimal.

    Return
    ------
    A string with comma-delimited values with the ID and mean and uncertainty
    for each item in the list. This is suitable for writing output to a .csv
    file.

    Example:
    --------
    import dtsa2.jmGen as jmg
    import gov.nist.microanalysis.Utility as epu
    nmZnO1   = 40.1
    uvOKa1   = epu.UncertainValue2(0.269157,0.000126)
    uvZnLa1  = epu.UncertainValue2(0.259251,9.4e-05)
    uvSiKa1  = epu.UncertainValue2(0.654561,8.4e-05)
    l_uvals = [uvOKa1, uvZnLa1, uvSiKa1]
    out = jmg.make_list_unc_val_string(nmZnO1, l_uvals)
    print(out)
    
    1> 40.1, 0.269157, 0.000126, 0.259251, 0.000094, 0.654561, 0.000084
    """
    lv = len(luv)
    i = 0
    rv = "%s, " % (id)
    for uv in luv:
        rc = round(uv.doubleValue(), 6)
        uc = round(uv.uncertainty(), 6)
        
        if i == lv-1:
            rv += "%g, %.6f" % (rc, uc)
        else:
            rv += "%g, %.6f, " % (rc, uc)
        i += 1
    return(rv)

print(help(jmg.make_list_unc_val_string))

l_uvals = [uvOKa1, uvZnLa1, uvSiKa1]

out = jmg.make_list_unc_val_string(nmZnO1, l_uvals)

print(out)




foo = epu.UncertainValue2(1.2345678, 0.0005678)
print(foo)

rc = foo.doubleValue()
uc = foo.uncertainty()

print(rc)
print(uc)

def pretty_print_unc_val(uv, n_digits):
    rc = round(uv.doubleValue(), n_digits)
    uc = round(uv.uncertainty(), n_digits)
    rv = u'%g \u00B1' % (rc)
    rv += ' %g' % (uc)
    return(rv)

print(pretty_print_unc_val(foo, 4))

# clean up cruft
shutil.rmtree(rptDir)
print "Done!"

end = time.time()
delta = (end-start)/60
msg = "This script required %.3f min" % delta
print msg
if(delta > 60):
    delta = delta/60
    msg = "...or %.3f hr" % delta
    print msg



