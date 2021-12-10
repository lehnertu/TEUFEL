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
parser.add_argument('file', help='the file name of the watch point HDF5 file')
parser.add_argument('freq', help="frequency in Hz", type=float)

print()
args = parser.parse_args()
radfile = args.file
radOK = os.path.isfile(radfile)
if not radOK:
  print("file not found")
  sys.exit()
f = args.freq

# Open the file for reading
print("reading ",radfile)
screen = TeufelScreen.read(radfile)
print()

(Nx,Ny) = screen.shape()
nots = screen.nots
dt = screen.dt

AmplitudeX = np.empty([Nx, Ny])
PhaseX = np.empty([Nx, Ny])
AmplitudeY = np.empty([Nx, Ny])
PhaseY = np.empty([Nx, Ny])
X = np.empty([Nx, Ny])
Y = np.empty([Nx, Ny])


for ix in range(Nx):
    for iy in range(Ny):
        # field data
        (X[ix,iy], Y[ix,iy]) = screen.screenPos(ix,iy)
        (EVec, BVec) = screen.getFieldTrace(ix,iy)
        t0 = screen.getStartTime(ix,iy)
        t = np.linspace(t0,t0+(nots-1)*dt,nots)
        cft = np.cos(2*np.pi*f*t)
        sft = np.sin(2*np.pi*f*t)
        # fourier components of the field
        Ex = EVec[:,0]
        cosEx = np.dot(Ex,cft)
        sinEx = np.dot(Ex,sft)
        cEx = cosEx + 1j*sinEx
        Ey = EVec[:,1]
        cosEy = np.dot(Ey,cft)
        sinEy = np.dot(Ey,sft)
        cEy = cosEy + 1j*sinEy
        Ez = EVec[:,2]
        cosEz = np.dot(Ez,cft)
        sinEz = np.dot(Ez,sft)
        cEz = cosEz + 1j*sinEz
        Bx = BVec[:,0]
        cosBx = np.dot(Bx,cft)
        sinBx = np.dot(Bx,sft)
        ccBx = cosBx - 1j*sinBx
        By = BVec[:,1]
        cosBy = np.dot(By,cft)
        sinBy = np.dot(By,sft)
        ccBy = cosBy - 1j*sinBy
        Bz = BVec[:,2]
        cosBz = np.dot(Bz,cft)
        sinBz = np.dot(Bz,sft)
        ccBz = cosBz - 1j*sinBz
        # Poynting vector S = (E.cross.B)/my0
        Sx = (cEy*ccBz - cEz*ccBy) *dt/mu0
        Sy = (cEz*ccBx - cEx*ccBz) *dt/mu0
        Sz = (cEx*ccBy - cEy*ccBx) *dt/mu0
        # Amplitude and Phase
        # Amplitude[ix,iy] = np.absolute(Sz)
        # Phase[ix,iy] = np.angle(Sz)
        AmplitudeX[ix,iy] = np.absolute(cEx)*dt
        PhaseX[ix,iy] = np.angle(cEx)
        AmplitudeY[ix,iy] = np.absolute(cEy)*dt
        PhaseY[ix,iy] = np.angle(cEy)



# figure with power density on screen

dX = screen.dx
dY = screen.dy
print("dx=%g dy=%g m" % (dX,dY))
# Etot = Amplitude.sum()*dX*dY
# print("integrated energy = ", 1e6*Etot, " µJ"(

maxX = np.max(AmplitudeX)
maxY = np.max(AmplitudeY)
maxV = np.max([maxX,maxY])
print("max = %g" % maxV)
expo = np.rint(np.log10(maxV))
scale = np.power(10,-expo)
print("scale = %g" % scale)

fig1 = plt.figure(1,figsize=(11,9))

ax1 = fig1.add_subplot(221)
plt.pcolormesh(X, Y, scale*AmplitudeX, cmap='CMRmap', vmin=0, vmax=scale*maxV)
plt.title('integrated $E_x$ amplitude [V/m s]')
plt.xlabel('x /m')
plt.ylabel('y /m')
cb=plt.colorbar()
cb.set_label(r'$10^{%d}$ Vs/m' % expo)

ax2 = fig1.add_subplot(222)
plt.pcolormesh(X, Y, PhaseX, cmap='seismic')
plt.title('$E_x$ phase [rad]')
plt.xlabel('x /m')
plt.ylabel('y /m')
cb=plt.colorbar()

ax3 = fig1.add_subplot(223)
plt.pcolormesh(X, Y, scale*AmplitudeY, cmap='CMRmap', vmin=0, vmax=scale*maxV)
plt.title('integrated $E_y$ amplitude [V/m s]')
plt.xlabel('x /m')
plt.ylabel('y /m')
cb=plt.colorbar()
cb.set_label(r'$10^{%d}$ Vs/m' % expo)

ax4 = fig1.add_subplot(224)
plt.pcolormesh(X, Y, PhaseY, cmap='seismic')
plt.title('$E_y$ phase [rad]')
plt.xlabel('x /m')
plt.ylabel('y /m')
cb=plt.colorbar()

fig1.tight_layout()

plt.show()
