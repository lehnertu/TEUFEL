#!/usr/bin/env python3

import sys, time
import os.path
import argparse
import numpy as np
import h5py
import sdds
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

pixels = 100

def PlotPS(x, y, xlabel='x', ylabel='y', rect_dens = [0.15, 0.1, 0.8, 0.8], center=True, square=False):
  global pixels
  # determine the scale range
  if center:
    if square:
        xmax=max(abs(x))
        ymax=max(abs(y))
        xmax=max([xmax,ymax])
        xmin=-xmax
        ymax=xmax
        ymin=-xmax
    else:
        xmax=max(abs(x))
        xmin=-xmax
        ymax=max(abs(y))
        ymin=-ymax
  else:
    xmin=min(x)
    xmax=max(x)
    ymin=min(y)
    ymax=max(y)
    if square:
        xmax=max([xmax,ymax])
        xmin=min([xmin,ymin])
        ymax=xmax
        ymin=xmin
  # Histogram in 2D
  if center:
    H, xticks, yticks = np.histogram2d(x,y,bins=pixels,
      range=[[-xmax,xmax],[-ymax,ymax]])
  else:
    H, xticks, yticks = np.histogram2d(x,y,bins=pixels,
      range=[[xmin,xmax],[ymin,ymax]])
  # wegen Matrixkonvention muss H transponiert werden
  # da von links oben beginnend gezeichnet wird, muss man vertikal flippen
  M = np.flipud(H.T)
  # erster Index wird von unten nach oben dargestellt
  # zweiter Index wird von links nach rechts dargestellt
  axDens = plt.axes(rect_dens)
  im = axDens.imshow(M, interpolation='nearest',
      aspect=(xticks[0]-xticks[-1])/(yticks[0]-yticks[-1]),
      extent=[xticks[0], xticks[-1], yticks[0], yticks[-1]])
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  return plt

def DistProp(X, prt=None):
    X = np.array(X)
    N = X.size
    avg = np.mean(X)
    dX = X - avg
    rms = np.sqrt(np.dot(dX,dX)/float(N))
    if prt != None:
        print(prt,"avg=%g  rms=%g" % (avg,rms))
    return (avg, rms)
    
parser = argparse.ArgumentParser()
parser.add_argument('file', help='the file name of the watch point file (HDF5 or SDDS)')
parser.add_argument('-pix', help="the number of pixels for the plots", dest="pix", type=int)
parser.add_argument('-img', help="output not to screen but image file", dest="img")

args = parser.parse_args()
if (args.pix != None): pixels = args.pix
bunfile = args.file
bunOK = os.path.isfile(bunfile)
if not bunOK:
  print("file not found")
  sys.exit()

# determine the file format
f = open(bunfile, "rb")
bytes = f.read(4)
is_sdds = (bytes == b'SDDS')
is_hdf5 = (bytes == b'\x89HDF')
if not (is_sdds or is_hdf5):
  print("file is neither SDDS nor HDF5 file")
  sys.exit()

# read HDF5 data
if is_hdf5:
    # Open the file for reading
    print("reading ",bunfile)
    hdf = h5py.File(bunfile, "r")
    print(hdf)
    print()

    # Get the group
    electrons = hdf['electrons']
    a = np.array(electrons)
    hdf.close()

    data = a.transpose()
    print(f'have read {data.shape} array.')

    x = data[0]
    y = data[1]
    z = data[2]
    bgx = data[3]
    bgy = data[4]
    bgz = data[5]

# read SDDS data
if is_sdds:
    # Open the file for reading
    print("reading ",bunfile)
    data = sdds.SDDS(0)
    data.load(bunfile)

    x = np.array(data.getColumnData("x"))
    xp = np.array(data.getColumnData("xp"))
    y = np.array(data.getColumnData("y"))
    yp = np.array(data.getColumnData("yp"))
    p = np.array(data.getColumnData("p"))
    z = np.array(data.getColumnData("z"))
    bgz = p*np.sqrt(1.0-np.square(xp)-np.square(yp))
    bgx = bgz*xp
    bgy = bgz*yp
    print(f'have read {x.shape} particles.')

dt = -(z-zmean) / 3.0e8
(t0, tau) = DistProp(dt, prt="arrival time t [s] : ")

p = np.sqrt(np.square(bgx)+np.square(bgy)+np.square(bgz))
(betagamma, p_rms) = DistProp(p, prt="particle momentum p [m_e*c] : ")
E = 0.511*np.sqrt(np.square(p)+1.0)
(E0, erms) = DistProp(E, prt="particle total energy E_tot [MeV] : ")
dE = E-E0
tau = 1.0e12 * tau
if tau!=0.0:
    chirp = 1.0e15 * ( np.dot(dE,dt)/float(Np) ) / (tau*tau)
else:
    chirp = 0.0
print()

dx = x-xmean
dy = y-ymean
xp = bgx/bgz
yp = bgy/bgz
DistProp(xp, prt="xp [] : ")
DistProp(yp, prt="yp [] : ")

dxp = (bgx-bgxmean)/bgz
dyp = (bgy-bgymean)/bgz
ex_rms = 1.0e6 * betagamma * np.sqrt((np.dot(dx,dx)*np.dot(dxp,dxp)-pow(np.dot(dx,dxp),2))/pow(float(Np),2))
ey_rms = 1.0e6 * betagamma * np.sqrt((np.dot(dy,dy)*np.dot(dyp,dyp)-pow(np.dot(dy,dyp),2))/pow(float(Np),2))
ez_rms = 1.0e15 * np.sqrt((np.dot(dE,dE)*np.dot(dt,dt)-pow(np.dot(dE,dt),2))/pow(float(Np),2))

string = r'$E_{beam}$ = %7.3f MeV' % E0 + '\n' + \
  '%d particles' % Np + '\n' + \
  r'$\sigma_E$ = %4.1f keV' % erms + '\n' + \
  r'$\sigma_t$ = %1.3f ps' % tau + '\n' + \
  r'$c$ = %2.2f keV/ps' % chirp + '\n' + \
  r'$\epsilon^n_{x, RMS}$ = %2.3f $\mu m$' % ex_rms  + '\n' + \
  r'$\epsilon^n_{y, RMS}$ = %2.3f $\mu m$' % ey_rms  + '\n' + \
  r'$\epsilon_{z, RMS}$ = %3.1f keV ps' % ez_rms
fig = plt.figure(figsize=(16,10),dpi=80)
# fig.text(0.78,0.9, string, color='b', size=20.0 , **text_params, linespacing=2)
fig.text(0.68,0.93, string, color='b', size=18.0 , ha='left', va='top', linespacing=2, family='serif')
PlotPS(1000.0*x,1000.0*xp,"$x$ [mm]", "$p_x$ [mrad]", rect_dens = [0.03, 0.55, 0.3, 0.4])
PlotPS(1000.0*y,1000.0*yp,"$y$ [mm]", "$p_y$ [mrad]", rect_dens = [0.35, 0.55, 0.3, 0.4])
PlotPS(1000.0*x,1000.0*y,"$x$ [mm]", "$y$ [mm]", rect_dens = [0.03, 0.08, 0.3, 0.4])
PlotPS(1e12*dt,E,"$dt$ [ps]", "$E [MeV]$", rect_dens = [0.35, 0.08, 0.3, 0.4], center=False)
PlotPS(1e12*dt,1000*dE-chirp*1e12*dt,"$dt$ [ps]", "$dE [keV]$", rect_dens = [0.68, 0.08, 0.3, 0.4], center=False)
# plt.annotate('E = ', xy=(0.8, 0.8),  xycoords='figure fraction',
#              xytext=(20, 20), textcoords='offset points',
#              ha="left", va="bottom",
#              )

if (args.img != None):
    plt.savefig(args.img)
else:
    plt.show()
