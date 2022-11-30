#!/usr/bin/env python3
# coding=UTF-8

import sys, time
import os.path
import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt
from screen import *
from matplotlib.ticker import NullFormatter
from matplotlib.patches import Circle

parser = argparse.ArgumentParser()
parser.add_argument('file', help='the file name of the HDF5 file of the screens recorded fields')
parser.add_argument('-fmax', help="plot range for the spectrum in Hz", dest="fmax", type=float)
parser.add_argument('-roi', help="ROI for the spectrum in Hz", dest="roi", type=float, nargs=2)
parser.add_argument('-circ', help="ROI for the spatial dist.", dest="circ", type=float)

print()
args = parser.parse_args()
if args.roi != None:
  roiOK = True
  f1 = args.roi[0]
  f2 = args.roi[1]
  print("frequency ROI : %g ... %g Hz" % (f1, f2))
else:
  roiOK = False

radfile = args.file
radOK = os.path.isfile(radfile)
if not radOK:
  print("file not found")
  sys.exit()

# Open the file for reading
print("reading ",radfile)
screen = TeufelScreen.read(radfile)
print()

print(f'peak Ex = {np.abs(screen.A[:,:,:,0]).max()} V/m')
print(f'peak Ey = {np.abs(screen.A[:,:,:,1]).max()} V/m')
print(f'peak Ez = {np.abs(screen.A[:,:,:,2]).max()} V/m')
print(f'peak Bx = {np.abs(screen.A[:,:,:,3]).max()} T')
print(f'peak By = {np.abs(screen.A[:,:,:,4]).max()} T')
print(f'peak Bz = {np.abs(screen.A[:,:,:,5]).max()} T')
print()

# analyze power density on axis
(EVec,BVec) = screen.getFieldTrace(screen.xcenter,screen.ycenter)
SVec = np.cross(EVec, BVec) / mu0
print('energy flow density on axis = ', 1e6*SVec.sum(axis=0)*screen.dt, " µJ/m²")
print()

# Figure 1 : spectral distribution
# Fig. 1a) spectrum on axis

df = 1.0/screen.dt/screen.nots
(f, amplit) = screen.spectralDensityH(screen.xcenter,screen.ycenter)
intspec = amplit.sum()*df
print("total spectral power density on axis = %g µJ/m²" % (1e6*intspec))

fig1 = plt.figure(1,figsize=(12,9))
ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8])
l1 = ax1.plot(1e-12*f, 1e6*1e12*amplit, "r-")
ax1.set_xlabel(r'$f$ [THz]')
ax1.set_ylabel(r'spectral power density on axis   $dI/dA/df$ [$\mu$J/(m$^2$ THz)]')
ax1.yaxis.label.set_color('r')
ax1.tick_params(axis='y', colors='r')
if roiOK:
  ax1.axvspan(1e-12*f1, 1e-12*f2, facecolor='#2ca02c', alpha=0.5)

if roiOK:
  nf1 = np.ceil(f1/df).astype('int')
  nf2 = np.floor(f2/df).astype('int')
  intspec = amplit[nf1:nf2].sum()*df
  print("integrated spectral power density in ROI = %g µJ/m²" % (1e6*intspec))

# Fig. 1b) spectrum integrated over the whole area
# and integrated over the cirular range (if defined)

totamp = np.zeros_like(amplit)
if args.circ != None:
    r2 = args.circ * args.circ
    circamp = np.zeros_like(amplit)
for ix in range(screen.Nx):
  for iy in range(screen.Ny):
    (f, amplit) = screen.spectralDensityH(ix,iy)  
    totamp += amplit
    if args.circ != None:
      (x,y) = screen.screenPos(ix,iy)
      if x*x + y*y <= r2:
        circamp += amplit
totamp *= screen.dx*screen.dy
if args.circ != None:
    circamp *= screen.dx*screen.dy

ax2 = ax1.twinx()
l2 = ax2.plot(1e-12*f, 1e6*1e12*totamp, "b-")
if args.circ != None:
    l3 = ax2.plot(1e-12*f, 1e6*1e12*circamp, "b--")
if args.fmax != None:
    ax2.set_xlim(0.0,1e-12*args.fmax)
ax2.set_ylabel(r'area integrated spectral power density   $dE/df$ [$\mu$J/THz]')
ax2.yaxis.label.set_color('b')
ax2.tick_params(axis='y', colors='b')

# Figure 2 : total power density distribution on screen

X = np.empty([screen.Nx, screen.Ny])
Y = np.empty([screen.Nx, screen.Ny])
Pz = np.empty([screen.Nx, screen.Ny])
for ix in range(screen.Nx):
  for iy in range(screen.Ny):
    (x,y) = screen.screenPos(ix,iy)
    X[ix,iy] = x
    Y[ix,iy] = y
    Pz[ix,iy] = screen.totalPowerDensity(ix,iy)
Etot = Pz.sum()*screen.dx*screen.dy
print("integrated energy = ", 1e6*Etot, " µJ")

if args.circ != None:
    Ecirc = 0
    r2 = args.circ * args.circ
    for ix in range(screen.Nx):
      for iy in range(screen.Ny):
        (x,y) = screen.screenPos(ix,iy)
        if x*x + y*y <= r2:
          Ecirc += Pz[ix,iy]
    Ecirc *= screen.dx*screen.dy
    print("integrated energy in circle = ", 1e6*Ecirc, " µJ")

fig2 = plt.figure(2,figsize=(12,9))
ax3 = fig2.add_subplot(111)
plt.contourf(X, Y, 1e6*Pz, 15, cmap='CMRmap')
plt.title('total energy density [$\mu$J/m$^2$]')
plt.xlabel('x /m')
plt.ylabel('y /m')
cb=plt.colorbar()
cb.set_label(r'energy density [$\mu$J/m$^2$]')
if args.circ != None:
    c = Circle(xy=[0,0], radius=args.circ)
    ax3.add_artist(c)
    c.set_fill(False)
    c.set_linestyle('dashed')
    c.set_lw(3)
    c.set_color('white')

# Figure 3 : power density distribution in ROI on screen

if roiOK:
    for ix in range(screen.Nx):
      for iy in range(screen.Ny):
        # we already have X and Y
        (f, amplit) = screen.spectralDensityH(ix,iy)
        Pz[ix,iy] = amplit[nf1:nf2].sum()*df
        
    Eroi = Pz.sum()*screen.dx*screen.dy
    print("integrated energy in ROI = ", 1e6*Eroi, " µJ")

    if args.circ != None:
        Ecirc = 0
        r2 = args.circ * args.circ
        for ix in range(screen.Nx):
          for iy in range(screen.Ny):
            (x,y) = screen.screenPos(ix,iy)
            if x*x + y*y <= r2:
              Ecirc += Pz[ix,iy]
        Ecirc *= screen.dx*screen.dy
        print("integrated energy in ROI in circle = ", 1e6*Ecirc, " µJ")

    fig3 = plt.figure(3,figsize=(12,9))
    ax4 = fig3.add_subplot(111)
    plt.contourf(X, Y, 1e6*Pz, 15, cmap='CMRmap')
    plt.title('energy density within ROI [$\mu$J/m$^2$]')
    plt.xlabel('x /m')
    plt.ylabel('y /m')
    cb=plt.colorbar()
    cb.set_label(r'energy density [$\mu$J/m$^2$]')
    if args.circ != None:
        c = Circle(xy=[0,0], radius=args.circ)
        ax4.add_artist(c)
        c.set_fill(False)
        c.set_linestyle('dashed')
        c.set_lw(3)
        c.set_color('white')

plt.show()

