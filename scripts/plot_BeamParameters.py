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

left, width = 0.15, 0.80
rect1 = [left, 0.55, width, 0.40]  #left, bottom, width, height
rect2 = [left, 0.08, width, 0.40]
fig = plt.figure(1,figsize=(12,9))

ax1 = fig.add_axes(rect1)
ax4 = fig.add_axes(rect2, sharex=ax1)

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
    
    l1 = ax1.plot(t, x, "r-", label=r'$x_{avg.}$')
    l2 = ax1.plot(t, y, "b-", label=r'$y_{avg.}$')
    # l3 = ax1.plot(t, z, "g-", label=r'$z_{avg.}$')

    ax1.set_ylabel(r'$position$ [m]')
    # lines = l1 + l2 + l3
    lines = l1 + l2
    labels = [l.get_label() for l in lines]
    ax1.legend(lines,labels,loc='upper right')
    for label in ax1.get_xticklabels():
	label.set_visible(False)
    ax1.grid(True)

    l4 = ax4.plot(t, bgx, "r-", label=r'$\beta_x\gamma$')
    l5 = ax4.plot(t, bgy, "b-", label=r'$\beta_y\gamma$')
    l6 = ax4.plot(t, bgz, "g-", label=r'$\beta_z\gamma$')
    
    ax4.set_ylabel(r'$\beta\gamma$')
    ax4.set_xlabel(r't [ns]')
    lines = l4 + l5 +l6
    labels = [l.get_label() for l in lines]
    ax4.legend(lines,labels,loc='upper right')
    ax4.grid(True)

plt.show()