#!/usr/bin/python

import sys, sdds, time
import os.path
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('file', help='the name of the SDDS file with a single-particle trajectory')
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

t = np.array(data.getColumnData("t"))*1e9
x = np.array(data.getColumnData("x"))*1e3
y = np.array(data.getColumnData("y"))*1e3
z = np.array(data.getColumnData("z"))*1e3
px = np.array(data.getColumnData("px"))
py = np.array(data.getColumnData("py"))
pz = np.array(data.getColumnData("pz"))

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

l1 = ax1.plot(t, x, "r-", label=r'x [mm]')
l2 = ax2.plot(t, y, "b-", label=r'y [mm]')
l3 = ax3.plot(t, z, "g-", label=r'z [mm]')
ax1.set_ylabel(r'x [mm]',color='Red')
ax1.tick_params(axis='y', colors='Red')
ax1.ticklabel_format(useOffset=False)
ax2.set_ylabel(r'y [mm]',color='Blue')
ax2.tick_params(axis='y', colors='Blue')
ax2.ticklabel_format(useOffset=False)
ax3.set_ylabel(r'z [mm]',color='Green')
ax3.tick_params(axis='y', colors='Green')
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

l4 = ax4.plot(t, px, "r-", label=r'$\beta\gamma_x$')
l5 = ax5.plot(t, py, "b-", label=r'$\beta\gamma_y$')
l6 = ax6.plot(t, pz, "g-", label=r'$\beta\gamma_z$')
ax4.set_ylabel(r'$\beta\gamma_x$',color='Red')
ax4.tick_params(axis='y', colors='Red')
ax4.ticklabel_format(useOffset=False)
ax5.set_ylabel(r'$\beta\gamma_y$',color='Blue')
ax5.tick_params(axis='y', colors='Blue')
ax5.ticklabel_format(useOffset=False)
ax6.set_ylabel(r'$\beta\gamma_z$',color='Green')
ax6.tick_params(axis='y', colors='Green')
ax6.ticklabel_format(useOffset=False)
ax4.set_xlabel(r't [ns]')
lines = l4 + l5 +l6
labels = [l.get_label() for l in lines]
ax6.legend(lines,labels,loc='upper right')
ax4.grid(True)

plt.show()
