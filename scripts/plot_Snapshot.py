#!/usr/bin/env python
# coding=UTF-8

import sys, time
import os.path
import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('file', help='the file name of the snapshot HDF5 file')
parser.add_argument('-Emax', help="plot range for the E field [V/m]", dest="emax", type=float)
parser.add_argument('-Bmax', help="plot range for the B field [T]", dest="bmax", type=float)

print
args = parser.parse_args()

radfile = args.file
radOK = os.path.isfile(radfile)
if not radOK:
  print "file not found"
  sys.exit()

# Open the file for reading
print "reading ",radfile
hdf = h5py.File(radfile, "r")
print hdf
print

# Get the groups
pos = hdf['ObservationPosition']
Nx = pos.attrs.get('Nx')
Ny = pos.attrs.get('Ny')
print "Nx=%d Ny=%d" % (Nx,Ny)
print pos
field = hdf['ElMagField']
print field
pos = np.array(pos)
a = np.array(field)
hdf.close()
print

xcenter = (Nx-1)/2
ycenter = (Ny-1)/2
centerposition = pos[xcenter][ycenter]
print "center position = ",centerposition

Ex = a[:,:,0]
Ey = a[:,:,1]
Ez = a[:,:,2]
Bx = a[:,:,3]
By = a[:,:,4]
Bz = a[:,:,5]

if args.emax == None:
    exmax = np.max(np.abs(Ex))
    eymax = np.max(np.abs(Ey))
    ezmax = np.max(np.abs(Ez))
    emax = np.max(np.array([exmax,eymax,ezmax]))
else:
    emax = args.emax
print "scale E = ",emax
elevels = np.linspace(-emax, emax, num=11)

if args.bmax == None:
    bxmax = np.max(np.abs(Bx))
    bymax = np.max(np.abs(By))
    bzmax = np.max(np.abs(Bz))
    bmax = np.max(np.array([bxmax,bymax,bzmax]))
else:
    bmax = args.bmax
print "scale B = ",bmax
blevels = np.linspace(-bmax, bmax, num=11)

fig1 = plt.figure(figsize=(12,8),dpi=80)
ax1 = plt.subplot(2,3,1)
ax2 = plt.subplot(2,3,2)
ax3 = plt.subplot(2,3,3)
ax4 = plt.subplot(2,3,4)
ax5 = plt.subplot(2,3,5)
ax6 = plt.subplot(2,3,6)
ax1.contourf(Ex,elevels,cmap='jet')
ax2.contourf(Ey,elevels,cmap='jet')
ax3.contourf(Ez,elevels,cmap='jet')
ax4.contourf(Bx,blevels,cmap='jet')
ax5.contourf(By,blevels,cmap='jet')
ax6.contourf(Bz,blevels,cmap='jet')

fig1.subplots_adjust(hspace=0.2, wspace=0.2)

plt.show()

"""
fig = plt.figure(1,figsize=(12,9))
ax3 = fig3.add_subplot(111)
plt.contourf(X, Y, 1e6*Pz, 15, cmap='CMRmap')
plt.title('total energy density [$\mu$J/m$^2$]')
plt.xlabel('x /m')
plt.ylabel('y /m')
cb=plt.colorbar()
cb.set_label(r'energy density [$\mu$J/m$^2$]')
"""
