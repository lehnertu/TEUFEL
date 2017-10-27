#!/usr/bin/python
# coding=UTF-8

import sys, time
import os.path
import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

# magnetic field constant in N/A²
mu0 = 4*np.pi*1e-7

parser = argparse.ArgumentParser()
parser.add_argument('file', help='the file name of the watch point HDF5 file')

args = parser.parse_args()
radfile = args.file
radOK = os.path.isfile(radfile)
if not radOK:
  print "file not found"
  sys.exit()

print "reading ",radfile

# Open the file for reading
hdf = h5py.File(radfile, "r")
print hdf
# Get the groups
pos = hdf['ObservationPosition']
Nx = pos.attrs.get('Nx')
Ny = pos.attrs.get('Ny')
print "Nx=%d Ny=%d" % (Nx,Ny)
print pos
field = hdf['ElMagField']
print field
t0 = field.attrs.get('t0')
dt = field.attrs.get('dt')
nots = field.attrs.get('NOTS')
print "t0=%g dt=%g NOTS=%d" % (t0, dt, nots)
pos = np.array(pos)
a = np.array(field)
hdf.close()

xcenter = (Nx-1)/2
ycenter = (Ny-1)/2
print "center = (",xcenter,",",ycenter,")"
centerposition = pos[xcenter][ycenter]
print "position = ",centerposition

onaxis = a[xcenter][ycenter]
data = onaxis.transpose()

Ex = data[0]
Ey = data[1]
Ez = data[2]
Bx = data[3]
By = data[4]
Bz = data[5]

EVec = np.array([Ex, Ey, Ez]).transpose()
BVec = np.array([Bx, By, Bz]).transpose()
# Poynting vector in V/m * (N/(A m)) / (N/A²) = W/m²
SVec = np.cross(EVec, BVec) / mu0
t = 1e9*np.arange(t0,t0+nots*dt,dt)
print 'energy flow density = ', SVec.sum(axis=0)*dt, " Ws/m²"

# first figure with the time-trace of the fields on axis

left, width = 0.15, 0.80
rect1 = [left, 0.55, width, 0.40]  #left, bottom, width, height
rect2 = [left, 0.08, width, 0.40]
fig = plt.figure(1,figsize=(12,9))

ax1 = fig.add_axes(rect1)
ax4 = fig.add_axes(rect2, sharex=ax1)

l1 = ax1.plot(t, Ex, "r-", label=r'$E_x$')
l2 = ax1.plot(t, Ey, "b-", label=r'$E_y$')
l3 = ax1.plot(t, Ez, "g-", label=r'$E_z$')

ax1.set_ylabel(r'$E$ [V/m]')
lines = l1 + l2 + l3
labels = [l.get_label() for l in lines]
ax1.legend(lines,labels,loc='upper right')
for label in ax1.get_xticklabels():
    label.set_visible(False)
ax1.grid(True)

l4 = ax4.plot(t, Bx, "r-", label=r'$B_x$')
l5 = ax4.plot(t, By, "b-", label=r'$B_y$')
l6 = ax4.plot(t, Bz, "g-", label=r'$B_z$')

ax4.set_ylabel(r'$B$ [T]')
ax4.set_xlabel(r't [ns]')
lines = l4 + l5 +l6
labels = [l.get_label() for l in lines]
ax4.legend(lines,labels,loc='upper right')
ax4.grid(True)

# second figure with power density on screen

t = 1e9*np.arange(t0,t0+nots*dt,dt)
X = np.empty([Nx, Ny])
Y = np.empty([Nx, Ny])
Pz = np.empty([Nx, Ny])
for ix in range(Nx):
    for iy in range(Ny):
	X[ix,iy] = pos[ix,iy,0]
	Y[ix,iy] = pos[ix,iy,1]
	trace = a[ix,iy]
	data = trace.transpose()
	Ex = data[0]
	Ey = data[1]
	Ez = data[2]
	Bx = data[3]
	By = data[4]
	Bz = data[5]
	EVec = np.array([Ex, Ey, Ez]).transpose()
	BVec = np.array([Bx, By, Bz]).transpose()
	SVec = np.cross(EVec, BVec) / mu0
	Pz[ix,iy] = (SVec.sum(axis=0))[2]*dt
print Pz

fig2 = plt.figure(2,figsize=(12,9))

plt.contourf(X, Y, Pz, 20, cmap='RdGy')
plt.colorbar();

plt.show()
