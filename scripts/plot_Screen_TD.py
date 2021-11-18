#!/usr/bin/env python3
# coding=UTF-8

import sys, time
import os.path
import argparse
import numpy as np
import h5py
from screen import *
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.patches import Circle

# magnetic field constant in N/A²
mu0 = 4*np.pi*1e-7

parser = argparse.ArgumentParser()
parser.add_argument('file', help='the file name of the screen output HDF5 file')
parser.add_argument('-xy', help="indeces of plot point", dest="xy", type=int, nargs=2)

print()
args = parser.parse_args()

radfile = args.file
radOK = os.path.isfile(radfile)
if not radOK:
  print("file not found")
  sys.exit()

# Open the file for reading
print("reading ",radfile)
screen = TeufelScreen.read(radfile)
print()

# analyze power density on axis
(EVec,BVec) = screen.getFieldTrace(screen.xcenter,screen.ycenter)
SVec = np.cross(EVec, BVec) / mu0
print('energy flow density on axis = ', 1e6*SVec.sum(axis=0)*screen.dt, " µJ/m²")
print()

# first figure with the time-trace of the fields on axis

t0 = screen.getStartTime(screen.xcenter,screen.ycenter)
t = 1e9*np.linspace(t0,t0+(screen.nots-1)*screen.dt,screen.nots)

left, width = 0.15, 0.80
rect1 = [left, 0.55, width, 0.40]  #left, bottom, width, height
rect2 = [left, 0.08, width, 0.40]
fig = plt.figure(1,figsize=(12,9))

ax1 = fig.add_axes(rect1)
ax2 = fig.add_axes(rect2, sharex=ax1)

l1 = ax1.plot(t, EVec[:,0], "r-", label=r'$E_x$')
l2 = ax1.plot(t, EVec[:,1], "b-", label=r'$E_y$')
l3 = ax1.plot(t, EVec[:,2], "g-", label=r'$E_z$')

ax1.set_ylabel(r'$E$ [V/m]')
lines = l1 + l2 + l3
labels = [l.get_label() for l in lines]
ax1.legend(lines,labels,loc='upper right')
for label in ax1.get_xticklabels():
    label.set_visible(False)
ax1.grid(True)

l4 = ax2.plot(t, BVec[:,0], "r-", label=r'$B_x$')
l5 = ax2.plot(t, BVec[:,1], "b-", label=r'$B_y$')
l6 = ax2.plot(t, BVec[:,2], "g-", label=r'$B_z$')

ax2.set_ylabel(r'$B$ [T]')
ax2.set_xlabel(r't [ns]')
lines = l4 + l5 +l6
labels = [l.get_label() for l in lines]
ax2.legend(lines,labels,loc='upper right')
ax2.grid(True)

# second figure with the time-trace of the fields off axis

if args.xy != None:

    xi = args.xy[0]
    yi = args.xy[1]
    print("index = (",xi,",",yi,")")
    position = screen.screenPos(xi,yi)
    print("off-axis position [m] = ",position)

    (EVec,BVec) = screen.getFieldTrace(xi,yi)
    SVec = np.cross(EVec, BVec) / mu0
    print('off axis energy flow density = ', 1e6*SVec.sum(axis=0)*screen.dt, " µJ/m²")
    
    t0 = screen.getStartTime(xi,yi)
    t = 1e9*np.linspace(t0,t0+(screen.nots-1)*screen.dt,screen.nots)

    fig2 = plt.figure(2,figsize=(12,9))
    ax3 = fig2.add_axes(rect1)
    ax4 = fig2.add_axes(rect2, sharex=ax3)
    
    l21 = ax3.plot(t, EVec[:,0], "r-", label=r'$E_x$')
    l22 = ax3.plot(t, EVec[:,1], "b-", label=r'$E_y$')
    l23 = ax3.plot(t, EVec[:,2], "g-", label=r'$E_z$')
    
    ax3.set_ylabel(r'$E$ [V/m]')
    lines = l21 + l22 + l23
    labels = [l.get_label() for l in lines]
    ax3.legend(lines,labels,loc='upper right')
    for label in ax3.get_xticklabels():
        label.set_visible(False)
    ax3.grid(True)

    l24 = ax4.plot(t, BVec[:,0], "r-", label=r'$B_x$')
    l25 = ax4.plot(t, BVec[:,1], "b-", label=r'$B_y$')
    l26 = ax4.plot(t, BVec[:,2], "g-", label=r'$B_z$')
    
    ax4.set_ylabel(r'$B$ [T]')
    ax4.set_xlabel(r't [ns]')
    lines = l24 + l25 +l26
    labels = [l.get_label() for l in lines]
    ax4.legend(lines,labels,loc='upper right')
    ax4.grid(True)

plt.show()
