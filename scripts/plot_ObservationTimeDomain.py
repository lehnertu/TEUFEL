#!/usr/bin/python

import sys, sdds, time
import os.path
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('files', default=[], nargs='+', help='file name(s) for the stored population(s)')
# parser.add_argument('file', help='the name of the SDDS file with the observation field data')
parser.add_argument('--list_params', dest='listpar',
  action='store_const', const=True, default=False, help='list all parameters available in the file')
parser.add_argument('--list_columns', dest='listcol',
  action='store_const', const=True, default=False, help='list all columns available in the file')

args = parser.parse_args()

left, width = 0.25, 0.70
rect1 = [left, 0.55, width, 0.40]  #left, bottom, width, height
rect2 = [left, 0.08, width, 0.40]
fig = plt.figure(1,figsize=(12,9))

ax1 = fig.add_axes(rect1)
ax2 = ax1.twinx()
ax3 = ax1.twinx()
ax4 = fig.add_axes(rect2, sharex=ax1)
ax5 = ax4.twinx()
ax6 = ax4.twinx()

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
    else:
	if ("t0" in data.parameterName and
	    "dt" in data.parameterName and
	    "NumberTimeSteps" in data.parameterName):
	    t0 = data.getParameterValue("t0")
	    dt = data.getParameterValue("dt")
	    nots = data.getParameterValue("NumberTimeSteps")
	    print r'%d steps starting at t0=%12.3g s step dt=%12.3g s' % (nots,t0,dt)
	    t = (np.linspace(t0, t0+nots*dt, num=nots)+0.5*dt)*1e9

    Ex = np.array(data.getColumnData("Ex"))
    Ey = np.array(data.getColumnData("Ey"))
    Ez = np.array(data.getColumnData("Ez"))
    Bx = np.array(data.getColumnData("Bx"))
    By = np.array(data.getColumnData("By"))
    Bz = np.array(data.getColumnData("Bz"))

    ax2.spines['left'].set_position(('axes', -0.12))
    ax2.spines["left"].set_visible(True)
    ax2.yaxis.set_label_position('left')
    ax2.yaxis.set_ticks_position('left')
    ax3.spines['left'].set_position(('axes', -0.24))
    ax3.spines["left"].set_visible(True)
    ax3.yaxis.set_label_position('left')
    ax3.yaxis.set_ticks_position('left')

    l1 = ax1.plot(t, Ex, "r-", label=r'$E_x$')
    l2 = ax2.plot(t, Ey, "b-", label=r'$E_y$')
    l3 = ax3.plot(t, Ez, "g-", label=r'$E_z$')

    ax1.set_ylabel(r'$E_x$ [V/m]',color='Red')
    ax1.tick_params(axis='Ey', colors='Red')
    ax1.ticklabel_format(useOffset=False)
    ax2.set_ylabel(r'$E_y$ [V/m]',color='Blue')
    ax2.tick_params(axis='Ey', colors='Blue')
    ax2.ticklabel_format(useOffset=False)
    ax3.set_ylabel(r'$E_z$ [V/m]',color='Green')
    ax3.tick_params(axis='Ez', colors='Green')
    ax3.ticklabel_format(useOffset=False)
    lines = l1 + l2 + l3
    labels = [l.get_label() for l in lines]
    ax3.legend(lines,labels,loc='upper right')
    for label in ax1.get_xticklabels():
	label.set_visible(False)
    ax1.grid(True)

    ax5.spines['left'].set_position(('axes', -0.12))
    ax5.spines["left"].set_visible(True)
    ax5.yaxis.set_label_position('left')
    ax5.yaxis.set_ticks_position('left')
    ax6.spines['left'].set_position(('axes', -0.24))
    ax6.spines["left"].set_visible(True)
    ax6.yaxis.set_label_position('left')
    ax6.yaxis.set_ticks_position('left')

    l4 = ax4.plot(t, Bx, "r-", label=r'$B_x$')
    l5 = ax5.plot(t, By, "b-", label=r'$B_y$')
    l6 = ax6.plot(t, Bz, "g-", label=r'$B_z$')
    ax4.set_ylabel(r'$B_x$ [T]',color='Red')
    ax4.tick_params(axis='Bx', colors='Red')
    ax4.ticklabel_format(useOffset=False)
    ax5.set_ylabel(r'$B_y$',color='Blue')
    ax5.tick_params(axis='By', colors='Blue')
    ax5.ticklabel_format(useOffset=False)
    ax6.set_ylabel(r'$B_z$',color='Green')
    ax6.tick_params(axis='z', colors='Green')
    ax6.ticklabel_format(useOffset=False)
    ax4.set_xlabel(r't [ns]')
    lines = l4 + l5 +l6
    labels = [l.get_label() for l in lines]
    ax6.legend(lines,labels,loc='upper right')
    ax4.grid(True)

plt.show()
