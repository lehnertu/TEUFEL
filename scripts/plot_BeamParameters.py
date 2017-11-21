#!/usr/bin/python
# coding=UTF-8

import sys, sdds, time
import os.path
import argparse
import numpy as np
import matplotlib.pyplot as plt

# magnetic field constant in N/A²
mu0 = 4*np.pi*1e-7

parser = argparse.ArgumentParser()
parser.add_argument('files', default=[], nargs='+', help='file name(s) for the stored population(s)')
# parser.add_argument('file', help='the name of the SDDS file with the observation field data')
parser.add_argument('--list_params', dest='listpar',
  action='store_const', const=True, default=False, help='list all parameters available in the file')
parser.add_argument('--list_columns', dest='listcol',
  action='store_const', const=True, default=False, help='list all columns available in the file')

args = parser.parse_args()

left, width = 0.10, 0.80
rect1 = [left, 0.68, width, 0.28]  #left, bottom, width, height
rect2 = [left, 0.37, width, 0.28]
rect3 = [left, 0.06, width, 0.28]

fig1 = plt.figure(1,figsize=(12,10))
ax1 = fig1.add_axes(rect1)
ax2 = ax1.twinx()	# separate axis for z
ax4 = fig1.add_axes(rect2, sharex=ax1)
ax5 = ax4.twinx()	# separate axis for z
ax21 = fig1.add_axes(rect3)
ax22 = ax21.twinx()	# separate axis for z

for file in args.files:

    fileOK = os.path.isfile(file)
    if not fileOK:
	print "file not found"
	sys.exit()

    print "reading ",file
    data = sdds.SDDS(0)
    data.load(file)
    if args.listpar:
	data.listParameters()
    if args.listcol:
	data.listColumns()

    if ("t" in data.columnName):
	t = np.array(data.getColumnData("t"))*1e9

    x = np.array(data.getColumnData("x_av"))
    y = np.array(data.getColumnData("y_av"))
    z = np.array(data.getColumnData("z_av"))
    bgx = np.array(data.getColumnData("bgx_av"))
    bgy = np.array(data.getColumnData("bgy_av"))
    bgz = np.array(data.getColumnData("bgz_av"))
    x_rms = np.array(data.getColumnData("x_rms"))
    y_rms = np.array(data.getColumnData("y_rms"))
    z_rms = np.array(data.getColumnData("z_rms"))
    
    l1 = ax1.plot(t, x, "r-", label=r'$x_{avg.}$')
    l2 = ax1.plot(t, y, "b-", label=r'$y_{avg.}$')
    l3 = ax2.plot(t, z, "g-", label=r'$z_{avg.}$')

    ax1.set_ylabel(r'position [m]')
    ax1.ticklabel_format(useOffset=False)
    for label in ax1.get_xticklabels():
	label.set_visible(False)
    ax2.set_ylabel(r'z position [m]', color="g")
    ax2.tick_params('y', colors='g')
    ax2.ticklabel_format(useOffset=False)
    for label in ax2.get_xticklabels():
	label.set_visible(False)
    lines = l1 + l2 + l3
    labels = [l.get_label() for l in lines]
    ax2.legend(lines,labels,loc='upper right')
    for label in ax1.get_xticklabels():
	label.set_visible(False)
    ax1.grid(True)

    l4 = ax4.plot(t, bgx, "r-", label=r'$\beta_x\gamma$')
    l5 = ax4.plot(t, bgy, "b-", label=r'$\beta_y\gamma$')
    l6 = ax5.plot(t, bgz, "g-", label=r'$\beta_z\gamma$')
    
    ax4.set_ylabel(r'$\beta\gamma$')
    ax4.ticklabel_format(useOffset=False)
    for label in ax4.get_xticklabels():
	label.set_visible(False)
    ax5.set_ylabel(r'$\beta_z\gamma$', color="g")
    ax5.tick_params('y', colors='g')
    ax5.ticklabel_format(useOffset=False)
    for label in ax5.get_xticklabels():
	label.set_visible(False)
    lines = l4 + l5 +l6
    labels = [l.get_label() for l in lines]
    ax5.legend(lines,labels,loc='upper right')
    ax4.grid(True)

    l21 = ax21.plot(t, 1e3*x_rms, "r-", label=r'$x_{r.m.s.}$')
    l22 = ax21.plot(t, 1e3*y_rms, "b-", label=r'$y_{r.m.s.}$')
    l23 = ax22.plot(t, 1e12*z_rms/3.0e8, "g-", label=r'$z_{r.m.s.}$')

    ax21.set_ylabel(r'beam size [mm]')
    ax21.set_xlabel(r't [ns]')
    ax21.ticklabel_format(useOffset=False)
    ax22.set_ylabel(r'bunch length [ps]', color="g")
    ax22.tick_params('y', colors='g')
    ax22.ticklabel_format(useOffset=False)
    lines = l21 + l22 +l23
    labels = [l.get_label() for l in lines]
    ax22.legend(lines,labels,loc='upper right')
    ax21.grid(True)


plt.show()
