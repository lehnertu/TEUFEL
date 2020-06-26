#!/usr/bin/env python3
# coding=UTF-8

import os.path
import numpy as np
from scipy import constants
from MeshedFields import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('file', help='the file name of the watch point HDF5 file')
parser.add_argument('-fmax', help="plot range for the spectrum in Hz", dest="fmax", type=float)

args = parser.parse_args()

radfile = args.file
radOK = os.path.isfile(radfile)
if not radOK:
  print("file not found")
  sys.exit()

screen = MeshedField.ReadMeshedField(radfile)

print()
print("%d points" % len(screen.points))
print("%d triangles" % len(screen.triangles))
area = screen.MeshArea()
print("total mesh area = %7.3f cm²" % (1.0e4*np.sum(area)))
normals = screen.MeshNormals()
average = np.sum(normals, axis=0)/screen.Np
print("screen normal = %s" % average)

area = screen.MeshArea()
S = [np.linalg.norm(screen.EnergyFlowVector(i)) for i in range(screen.Np)]
Pz = [screen.NormalEnergyFlow(i) for i in range(screen.Np)]

print("peak energy density = %.6g J/m²" % np.max(S))
print("total pulse energy = %.6g µJ" % (1e6*np.dot(area,Pz)))

# plot figure with spectrum integrated over the whole screen

nf = screen.Nt
fmax = 1.0/screen.dt
f = np.linspace(0.0,fmax,nf)[:nf//2]
df=1.0/screen.dt/nf

print()
print("trace length = %d" % screen.Nt)
print("trace step size = %g s" % screen.dt)
print("FFT max. frequency = %g Hz" % fmax)
print("FFT bin width = %g Hz" % df)
print()

amplit = np.zeros(nf//2)
for index in range(screen.Np):
    trace = screen.FieldTrace(index)
    Ex = trace[:,0]
    Ey = trace[:,1]
    Bx = trace[:,3]
    By = trace[:,4]
    spectEx = np.fft.fft(Ex)[:nf//2]
    spectEy = np.fft.fft(Ey)[:nf//2]
    spectBx = np.fft.fft(Bx)[:nf//2]
    spectBy = np.fft.fft(By)[:nf//2]
    amplit += ( np.real(spectEx*np.conj(spectBy)/constants.mu_0) -
                np.real(spectEy*np.conj(spectBx)/constants.mu_0) ) * 2*screen.dt/(df*nf) * area[index]
# intspec = amplit.sum()*df
intspec = np.trapz(amplit)*df
print("integrated spectrum = %g µJ" % (1e6*intspec))

fig1 = plt.figure(1,figsize=(8,6))
ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8])
l1 = ax1.plot(1e-12*f, 1e6*1e12*amplit, "r-")
ax1.set_xlabel(r'$f$ [THz]')
ax1.set_ylabel(r'area integrated spectral power density   $dE/df$ [$\mu$J/THz]')
ax1.yaxis.label.set_color('r')
ax1.tick_params(axis='y', colors='r')
ax1.grid(True)
if args.fmax != None:
    ax1.set_xlim(0.0,1e-12*args.fmax)

"""
if roiOK:
  ax1.axvspan(1e-12*f1, 1e-12*f2, facecolor='#2ca02c', alpha=0.5)
"""

plt.show()


"""
if roiOK:
  nf1 = np.ceil(f1/df).astype('int')
  nf2 = np.floor(f2/df).astype('int')
  intspec = amplit[nf1:nf2].sum()*df
  print "integrated spectral power density in ROI = %g µJ/m²" % (1e6*intspec)

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
"""

# generate a VTK visualization of the power density

print()
print("VTK display")
print("-----------")
print("rotate scene with left mouse button + drag")
print("shift view with middle mouse button + drag")
print("zoom with mouse wheel")
print("click on cell to plot waveform")
print("close waveform window before next interaction")
print()

def pick(id):
    if id>0 and id<screen.Np:
        print("cell No. %d pos=%s" % (id,screen.pos[id]))
        print("pointing vector S=%s" % screen.EnergyFlowVector(id))
        screen.ShowFieldTrace(id)

screen.ShowMeshedField(scalars=Pz,scalarTitle="Pz",pickAction=pick,showGrid=False)


