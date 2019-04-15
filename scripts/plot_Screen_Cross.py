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
parser.add_argument('-x', help="position index for plotting a vertical cross section", dest="xindex", type=int)
parser.add_argument('-y', help="position index for plotting a horizontal cross section", dest="yindex", type=int)

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

if args.yindex != None:
    ply = args.yindex
    plotdir = 'x'
    print "plotting at iy=%d y=%f m" % (ply,pos[0,ply,1])
else:
    if args.xindex != None:
        plx = args.xindex
        plotdir = 'y'
        print "plotting at ix=%d x=%f m" % (plx,pos[plx,0,0])
    else:
        ply = ycenter
        plotdir = 'x'
        print "plotting at iy=%d y=%f m" % (ply,pos[0,ply,1])

# first figure with spectrum on axis

onaxis = a[xcenter][ycenter]
data = onaxis.transpose()

Ex = data[0]
Ey = data[1]
Ez = data[2]
Bx = data[3]
By = data[4]
Bz = data[5]

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
#print "total power density on axis = %g µJ/m²" % (1e6*intspec)

if roiOK:
  nf1 = np.ceil(f1/df).astype('int')
  nf2 = np.floor(f2/df).astype('int')
  intspec = amplit[nf1:nf2].sum()*df
  print "integrated spectral power density in ROI = %g µJ/m²" % (1e6*intspec)

# add spectrum integrated over the whole area

totamp = np.zeros_like(amplit)
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
    spectE = np.fft.fft(Ex)[:nots/2]
    spectB = np.fft.fft(By)[:nots/2]
    totamp += np.real(spectE*np.conj(spectB)/mu0)*2*dt/(df*nf)
dX = pos[1,0,0]-pos[0,0,0]
dY = pos[0,1,1]-pos[0,0,1]
totamp *= dX*dY
print

ax2 = ax1.twinx()
l2 = ax2.plot(1e-12*f, 1e6*1e12*totamp, "b-")
if args.fmax != None:
    ax2.set_xlim(0.0,1e-12*args.fmax)
ax2.set_ylabel(r'area integrated spectral power density   $dE/df$ [$\mu$J/THz]')
ax2.yaxis.label.set_color('b')
ax2.tick_params(axis='y', colors='b')

# plot a cross section
# collect spectra for every point

if plotdir == 'y':
    Np = Ny
    xp = [pos[plx,i,1] for i in range(Ny)]
    dp = dY
else:
    Np = Nx
    xp = [pos[i,ply,0] for i in range(Nx)]
    dp = dX
totamp = np.zeros(Np)
roiamp = np.zeros(Np)
Nf = np.ceil(nots/2).astype('int')
fp = np.arange(0.0,df*Nf,df)
if args.fmax != None:
    Nf = np.ceil(args.fmax/df).astype('int')
    fp = fp[0:Nf]
specdens = np.zeros((Np,Nf))
for ip in range(Np):
    if plotdir == 'y':
        trace = a[plx,ip]
    else:
        trace = a[ip,ply]
    data = trace.transpose()
    Ex = data[0]
    Ey = data[1]
    Ez = data[2]
    Bx = data[3]
    By = data[4]
    Bz = data[5]
    spectE = np.fft.fft(Ex)[:nots/2]
    spectB = np.fft.fft(By)[:nots/2]
    amplit = np.real(spectE*np.conj(spectB)/mu0)*2*dt/(df*nf)
    specdens[ip] = np.array(amplit)[0:Nf]
    totamp[ip] = amplit.sum()*df
    if roiOK:
        nf1 = np.ceil(f1/df).astype('int')
        nf2 = np.floor(f2/df).astype('int')
        roiamp[ip] = amplit[nf1:nf2].sum()*df

# gaussian fit

def gauss(x,mean,rms):
    return np.exp(-np.square((x-mean)/xrms)/2.0)/(np.sqrt(2.0*np.pi)*rms)
    
if roiOK:
    print "gaussian fit:"
    si = roiamp.sum()
    six = (roiamp*xp).sum()
    xcenter = six/si
    print "center = %f m" % xcenter
    dx = xp-xcenter
    six2 = (roiamp*dx*dx).sum()
    xrms = np.sqrt(six2/si)
    print "rms = %f m" % xrms
    
fig2 = plt.figure(2,figsize=(12,9))
ax1 = fig2.add_axes([0.1, 0.1, 0.8, 0.8])
l1 = ax1.plot(xp, 1e6*totamp, "b-")
if roiOK:
    l2 = ax1.plot(xp, 1e6*roiamp, "r-")
    xg = np.linspace(xp[0],xp[-1],100)
    fg = gauss(xg,xcenter,xrms)
    l3 = ax1.plot(xg, 1e6*si*dp*fg, "r--")
    
ax1.set_xlabel(r'position [m]')
ax1.set_ylabel(r'power density $dW/dA$ [$\mu$J/m$^2$]')

# transversal spectral density plot

fig3 = plt.figure(3,figsize=(12,9))
ax3 = fig3.add_subplot(111)
plt.contourf(fp, xp, 1e6*1e12*specdens, 15, cmap='CMRmap')
plt.title('spectral energy density [$\mu$J/(m$^2$ THz]')
plt.xlabel('f [THz]')
plt.ylabel('pos [m]')
cb=plt.colorbar()
cb.set_label(r'spectral energy density [$\mu$J/(m$^2$ THz]')

plt.show()
