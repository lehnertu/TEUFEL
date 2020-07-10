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
parser.add_argument('-roi', help="ROI for the spectrum in Hz", dest="roi", type=float, nargs=2)

args = parser.parse_args()

radfile = args.file
radOK = os.path.isfile(radfile)
if not radOK:
  print("file not found")
  sys.exit()

# if there is a ROI given we will plot the intensity in the ROI
# instead of the total intensity
if args.roi != None:
  roiOK = True
  f1 = args.roi[0]
  f2 = args.roi[1]
  print("frequency ROI : %g ... %g Hz" % (f1, f2))
else:
  roiOK = False

# --------------------------------------------------
# read and analyze the file
# --------------------------------------------------

screen = MeshedField.ReadMeshedField(radfile)

print()
print("%d points" % len(screen.points))
print("%d triangles" % len(screen.triangles))
area = screen.MeshArea()
print("total mesh area = %7.3f cm²" % (1.0e4*np.sum(area)))
normals = screen.MeshNormals()
average_normal = np.sum(normals, axis=0)/screen.Np
print("screen normal = %s" % average_normal)

area = screen.MeshArea()
S = [np.linalg.norm(screen.EnergyFlowVector(i)) for i in range(screen.Np)]
Pz = [screen.NormalEnergyFlow(i) for i in range(screen.Np)]

print("peak energy density = %.6g J/m²" % np.max(S))
print("total pulse energy = %.6g µJ" % (1e6*np.dot(area,Pz)))

# --------------------------------------------------
# compute and integrate the spectrum over the screen
# --------------------------------------------------

nf = screen.Nt
fmax = 1.0/screen.dt
f = np.linspace(0.0,fmax,nf)[:nf//2]
df=1.0/screen.dt/nf
print()
print("trace length = %d" % screen.Nt)
print("trace step size = %g s" % screen.dt)
print("FFT max. frequency = %g Hz" % fmax)
print("FFT bin width = %g Hz" % df)
if roiOK:
    nf1 = np.ceil(f1/df).astype('int')
    nf2 = np.floor(f2/df).astype('int')
    print("FFT ROI = [%d ... %d]" % (nf1,nf2))
print()

amp_tot = np.zeros(nf//2)
if roiOK: Pz_roi = np.zeros_like(Pz)

# loop over the triangular cells of the screen
for index in range(screen.Np):
    trace = screen.FieldTrace(index)
    Ex = trace[:,0]
    Ey = trace[:,1]
    Ez = trace[:,2]
    Bx = trace[:,3]
    By = trace[:,4]
    Bz = trace[:,5]
    # compute the spectra of the field components
    spectEx = np.fft.fft(Ex)[:nf//2]
    spectEy = np.fft.fft(Ey)[:nf//2]
    spectEz = np.fft.fft(Ez)[:nf//2]
    spectBx = np.fft.fft(Bx)[:nf//2]
    spectBy = np.fft.fft(By)[:nf//2]
    spectBz = np.fft.fft(Bz)[:nf//2]
    # compute the spectral components of the Poyntig vector
    amp_x = np.real(spectEy*np.conj(spectBz)/constants.mu_0) - np.real(spectEz*np.conj(spectBy)/constants.mu_0)
    amp_y = np.real(spectEz*np.conj(spectBx)/constants.mu_0) - np.real(spectEx*np.conj(spectBz)/constants.mu_0)
    amp_z = np.real(spectEx*np.conj(spectBy)/constants.mu_0) - np.real(spectEy*np.conj(spectBx)/constants.mu_0)
    # compute the component of the Poyntig vector against the direction of the average normal
    amp_norm = np.array([ np.dot( np.array([ax,ay,az]), -average_normal ) for (ax,ay,az) in zip(amp_x,amp_y,amp_z)])
    amp_tot += amp_norm * 2*screen.dt/(df*nf) * area[index]
    Pz_roi[index] = amp_norm[nf1:nf2].sum() * 2*screen.dt/(df*nf) * area[index]

print("integrated spectrum = %g µJ" % (1e6*np.trapz(amp_tot)*df))
if roiOK: print("integrated spectrum in ROI = %g µJ" % (1e6*Pz_roi.sum()*df))

# --------------------------------------------------
# plot the spectrum
# --------------------------------------------------

fig1 = plt.figure(1,figsize=(8,6))
ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8])
l1 = ax1.plot(1e-12*f, 1e6*1e12*amp_tot, "r-")
ax1.set_xlabel(r'$f$ [THz]')
ax1.set_ylabel(r'area integrated spectral power density   $dE/df$ [$\mu$J/THz]')
ax1.yaxis.label.set_color('r')
ax1.tick_params(axis='y', colors='r')
ax1.grid(True)
if args.fmax != None:
    ax1.set_xlim(0.0,1e-12*args.fmax)
if roiOK:
  ax1.axvspan(1e-12*f1, 1e-12*f2, facecolor='#2ca02c', alpha=0.5)
plt.show()

# --------------------------------------------------
# generate a VTK visualization of the power density
# --------------------------------------------------

print()
print("VTK display")
print("-----------")
print("rotate scene with left mouse button + drag")
print("shift view with middle mouse button + drag")
print("zoom with mouse wheel")
print("click on cell to plot waveform")
print("close waveform window before next click")
print()

def pick(id):
    if id>0 and id<screen.Np:
        print("cell No. %d pos=%s" % (id,screen.pos[id]))
        print("pointing vector S=%s" % screen.EnergyFlowVector(id))
        screen.ShowFieldTrace(id)

if roiOK:
    screen.ShowMeshedField(scalars=Pz_roi,scalarTitle="Pz(ROI)",pickAction=pick,showGrid=False)
else:
    screen.ShowMeshedField(scalars=Pz,scalarTitle="Pz",pickAction=pick,showGrid=False)


