#!/usr/bin/env python

import sys, sdds, time
import os.path
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('file', help='the name of the SDDS file with the observation field data')
parser.add_argument('--list_columns', dest='listcol',
  action='store_const', const=True, default=False, help='list all columns available in the file')

args = parser.parse_args()

fileOK = os.path.isfile(args.file)
if not fileOK:
  print "file not found"
  sys.exit()

print "reading ",args.file
data = sdds.SDDS(0)
data.load(args.file)
if args.listcol:
  data.listColumns()

f = np.array(data.getColumnData("f"))*1.0e-12
Ax = np.array(data.getColumnData("Ax"))
Ay = np.array(data.getColumnData("Ay"))
Az = np.array(data.getColumnData("Az"))
Px = np.array(data.getColumnData("Px"))
Py = np.array(data.getColumnData("Py"))
Pz = np.array(data.getColumnData("Pz"))

left, width = 0.25, 0.70
rect1 = [left, 0.55, width, 0.40]  #left, bottom, width, height
rect2 = [left, 0.08, width, 0.40]
fig = plt.figure(1,figsize=(12,9))

ax1 = fig.add_axes(rect1)
ax2 = ax1.twinx()
ax3 = ax1.twinx()

ax2.spines['left'].set_position(('axes', -0.12))
ax2.spines["left"].set_visible(True)
ax2.yaxis.set_label_position('left')
ax2.yaxis.set_ticks_position('left')
ax3.spines['left'].set_position(('axes', -0.24))
ax3.spines["left"].set_visible(True)
ax3.yaxis.set_label_position('left')
ax3.yaxis.set_ticks_position('left')

l1 = ax1.plot(f, Ax, "r-", label=r'$A_x$')
l2 = ax2.plot(f, Ay, "b-", label=r'$A_y$')
l3 = ax3.plot(f, Az, "g-", label=r'$A_z$')
ax1.set_ylabel(r'$A_x$',color='Red')
ax1.tick_params(axis='Ax', colors='Red')
ax1.ticklabel_format(useOffset=False)
ax2.set_ylabel(r'$A_y$',color='Blue')
ax2.tick_params(axis='Ay', colors='Blue')
ax2.ticklabel_format(useOffset=False)
ax3.set_ylabel(r'$A_z$',color='Green')
ax3.tick_params(axis='Az', colors='Green')
ax3.ticklabel_format(useOffset=False)
lines = l1 + l2 +l3
labels = [l.get_label() for l in lines]
ax3.legend(lines,labels,loc='upper right')
for label in ax1.get_xticklabels():
  label.set_visible(False)
ax1.grid(True)

ax4 = fig.add_axes(rect2, sharex=ax1)
ax5 = ax4.twinx()
ax6 = ax4.twinx()

ax5.spines['left'].set_position(('axes', -0.12))
ax5.spines["left"].set_visible(True)
ax5.yaxis.set_label_position('left')
ax5.yaxis.set_ticks_position('left')
ax6.spines['left'].set_position(('axes', -0.24))
ax6.spines["left"].set_visible(True)
ax6.yaxis.set_label_position('left')
ax6.yaxis.set_ticks_position('left')

l4 = ax4.plot(f, Px, "r-", label=r'$P_x$')
l5 = ax5.plot(f, Py, "b-", label=r'$P_y$')
l6 = ax6.plot(f, Pz, "g-", label=r'$P_z$')
ax4.set_ylabel(r'$P_x$ [rad]',color='Red')
ax4.tick_params(axis='Px', colors='Red')
ax4.ticklabel_format(useOffset=False)
ax5.set_ylabel(r'$P_y$ [rad]',color='Blue')
ax5.tick_params(axis='Py', colors='Blue')
ax5.ticklabel_format(useOffset=False)
ax6.set_ylabel(r'$P_z [rad]$',color='Green')
ax6.tick_params(axis='Pz', colors='Green')
ax6.ticklabel_format(useOffset=False)
ax4.set_xlabel(r'f [THz]')
lines = l4 + l5 +l6
labels = [l.get_label() for l in lines]
ax6.legend(lines,labels,loc='upper right')
ax4.grid(True)

plt.show()
