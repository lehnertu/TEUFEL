#!/usr/bin/env python
# coding=UTF-8

import sys, time
import os.path
import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.patches import Circle

# magnetic field constant in N/A²
mu0 = 4*np.pi*1e-7

parser = argparse.ArgumentParser()
parser.add_argument('file', help='the file name of the watch point HDF5 file')
parser.add_argument('-fmax', help="plot range for the spectrum in Hz", dest="fmax", type=float)
parser.add_argument('-roi', help="ROI for the spectrum in Hz", dest="roi", type=float, nargs=2)
parser.add_argument('-circ', help="ROI for the spatial dist.", dest="circ", type=float)

print
args = parser.parse_args()
if args.roi != None:
  roiOK = True
  f1 = args.roi[0]
  f2 = args.roi[1]
  print "frequency ROI : %g ... %g Hz" % (f1, f2)
else:
  roiOK = False

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
t0 = field.attrs.get('t0')
dt = field.attrs.get('dt')
nots = field.attrs.get('NOTS')
print "t0=%g dt=%g NOTS=%d" % (t0, dt, nots)
pos = np.array(pos)
a = np.array(field)
hdf.close()
print

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
# t = 1e9*np.arange(t0,t0+(nots-1)*dt,dt)
t = 1e9*np.linspace(t0,t0+(nots-1)*dt,nots)
print 'energy flow density = ', 1e6*SVec.sum(axis=0)*dt, " µJ/m²"

# second figure with spectrum on axis

nf = nots
fmax = 1.0/dt
f = np.linspace(0.0,fmax,nf)[:nots/2]
df=1.0/dt/nf
spectE = np.fft.fft(Ex)[:nots/2]
# print spectE[:30]
spectB = np.fft.fft(By)[:nots/2]
# print spectB[:30]
amplit = np.real(spectE*np.conj(spectB)/mu0)*2*dt/(df*nf)

fig1 = plt.figure(1,figsize=(12,9))
ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8])
l1 = ax1.plot(1e-12*f, 1e6*1e12*amplit, "r-")
ax1.set_xlabel(r'$f$ [THz]')
ax1.set_ylabel(r'spectral power density on axis   $dI/dA/df$ [$\mu$J/(m$^2$ THz)]')
ax1.yaxis.label.set_color('r')
ax1.tick_params(axis='y', colors='r')
if roiOK:
  ax1.axvspan(1e-12*f1, 1e-12*f2, facecolor='#2ca02c', alpha=0.5)

print
print "trace length = %d" % nots
print "trace step size = %g s" % dt
print "FFT max. frequency = %g Hz" % fmax
print "FFT bin width = %g Hz" % df
print
intspec = amplit.sum()*df
print "total spectral power density on axis = %g µJ/m²" % (1e6*intspec)

if roiOK:
  nf1 = np.ceil(f1/df).astype('int')
  nf2 = np.floor(f2/df).astype('int')
  intspec = amplit[nf1:nf2].sum()*df
  print "integrated spectral power density in ROI = %g µJ/m²" % (1e6*intspec)

# add spectrum integrated over the whole area
# and integrated over the cirular range (if defined)
totamp = np.zeros_like(amplit)
if args.circ != None:
    r2 = args.circ * args.circ
    circamp = np.zeros_like(amplit)
for ix in range(Nx):
  for iy in range(Ny):
    trace = a[ix,iy]
    data = trace.transpose()
    Ex = data[0]
    Ey = data[1]
    Ez = data[2]
    Bx = data[3]
    By = data[4]
    Bz = data[5]
    spectEx = np.fft.fft(Ex)[:nots/2]
    spectEy = np.fft.fft(Ey)[:nots/2]
    spectEz = np.fft.fft(Ez)[:nots/2]
    spectBx = np.fft.fft(Bx)[:nots/2]
    spectBy = np.fft.fft(By)[:nots/2]
    spectBz = np.fft.fft(Bz)[:nots/2]
    totamp += np.real(spectEx*np.conj(spectBy)/mu0)*2*dt/(df*nf)
    if args.circ != None:
      x = pos[ix,iy,0]
      y = pos[ix,iy,1]
      if x*x + y*y <= r2:
        circamp += np.real(spectEx*np.conj(spectBy)/mu0)*2*dt/(df*nf)

dX = pos[1,0,0]-pos[0,0,0]
dY = pos[0,1,1]-pos[0,0,1]
print "dx=%g dy=%g m" % (dX,dY)
totamp *= dX*dY
if args.circ != None:
    circamp *= dX*dY

ax2 = ax1.twinx()
l2 = ax2.plot(1e-12*f, 1e6*1e12*totamp, "b-")
if args.circ != None:
    l3 = ax2.plot(1e-12*f, 1e6*1e12*circamp, "b--")
if args.fmax != None:
    ax2.set_xlim(0.0,1e-12*args.fmax)
ax2.set_ylabel(r'area integrated spectral power density   $dE/df$ [$\mu$J/THz]')
ax2.yaxis.label.set_color('b')
ax2.tick_params(axis='y', colors='b')

# third figure with power density on screen

t = 1e9*np.linspace(t0,t0+(nots-1)*dt,nots)
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

dX = pos[1,0,0]-pos[0,0,0]
dY = pos[0,1,1]-pos[0,0,1]
print "dx=%g dy=%g m" % (dX,dY)
Etot = Pz.sum()*dX*dY
print "integrated energy = ", 1e6*Etot, " µJ"

if args.circ != None:
    Ecirc = 0
    r2 = args.circ * args.circ
    for ix in range(Nx):
      for iy in range(Ny):
	    x = pos[ix,iy,0]
	    y = pos[ix,iy,1]
	    if x*x + y*y <= r2:
	      Ecirc += Pz[ix,iy]
    Ecirc *= dX*dY
    print "integrated energy in circle = ", 1e6*Ecirc, " µJ"

fig3 = plt.figure(3,figsize=(12,9))
ax3 = fig3.add_subplot(111)
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


# fourth figure with power density in ROI on screen

if roiOK:
    t = 1e9*np.arange(t0,t0+(nots-1)*dt,dt)
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
        spectEx = np.fft.fft(Ex)[:nots/2]
        spectEy = np.fft.fft(Ey)[:nots/2]
        spectEz = np.fft.fft(Ez)[:nots/2]
        spectBx = np.fft.fft(Bx)[:nots/2]
        spectBy = np.fft.fft(By)[:nots/2]
        spectBz = np.fft.fft(Bz)[:nots/2]
        amplit = np.real(spectEx*np.conj(spectBy)-spectEy*np.conj(spectBx))/mu0*2*dt/(df*nf)
        Pz[ix,iy] = amplit[nf1:nf2].sum()*df

    dX = pos[1,0,0]-pos[0,0,0]
    dY = pos[0,1,1]-pos[0,0,1]
    Eroi = Pz.sum()*dX*dY
    print "integrated energy in ROI = ", 1e6*Eroi, " µJ"

    if args.circ != None:
        Ecirc = 0
        r2 = args.circ * args.circ
        for ix in range(Nx):
          for iy in range(Ny):
            x = pos[ix,iy,0]
            y = pos[ix,iy,1]
            if x*x + y*y <= r2:
              Ecirc += Pz[ix,iy]
        Ecirc *= dX*dY
        print "integrated energy in ROI in circle = ", 1e6*Ecirc, " µJ"

    fig4 = plt.figure(4,figsize=(12,9))
    ax4 = fig4.add_subplot(111)
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
