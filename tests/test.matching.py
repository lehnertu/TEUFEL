#!/usr/bin/env python3
# coding=UTF-8

import os, sys
import numpy as np

TEUFEL_ROOT = os.getcwd()
sys.path.append(TEUFEL_ROOT+'/scripts')

import sdds

data = sdds.SDDS(0)
data.load("matching.sdds")

z_av = np.array(data.getColumnData("z_av"))
y_rms = np.array(data.getColumnData("y_rms"))
y_test = np.array([y for y,z in zip(y_rms,z_av) if (z>0.5 and z<1.5)])

# known a priori
y_matched = 0.2359e-3
min_y = np.min(y_test)
max_y = np.max(y_test)

print("minimum = %6.3f mm" % (1e3*min_y) )
print("known   = %6.3f mm" % (1e3*y_matched) )
print("maximum = %6.3f mm" % (1e3*max_y) )

errors = 0

if (y_matched-min_y) > 0.05*y_matched : errors += 1
if (max_y-y_matched) > 0.05*y_matched : errors += 1
if (max_y-min_y) > 0.1*y_matched : errors += 1

if errors==0:
    print("\033[1;32mOK.\033[0m")
else:
    print("\033[1;31mERROR.\033[0m")

sys.exit(errors)
