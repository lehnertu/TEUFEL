#!/usr/bin/env python3

import sys, sdds, time
import os.path
import argparse
import numpy as np
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
  # erster Index ist jetzt Energie, wird von unten nach oben dargestellt
  # zweiter Index ist Zeit, wird von links nach rechts dargestellt
  axDens = plt.axes(rect_dens)
  im = axDens.imshow(M, interpolation='nearest',
      aspect=(xticks[0]-xticks[-1])/(yticks[0]-yticks[-1]),
      extent=[xticks[0], xticks[-1], yticks[0], yticks[-1]])
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  # plt.colorbar()
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
parser.add_argument('file', help='the file name of the watch point')
parser.add_argument('--list_columns', dest='listcol',
  action='store_const', const=True, default=False,
  help='list all columns available in the output files')
parser.add_argument('--list_params', dest='listpar',
  action='store_const', const=True, default=False, help='list all parameters available in the file')
parser.add_argument('--cp', dest='currpro',
  action='store_const', const=True, default=False,
  help='plot a current profile')
parser.add_argument('-pix', help="the number of pixels for the plots", dest="pix", type=int)
parser.add_argument('-img', help="output not to screen but image file", dest="img")

args = parser.parse_args()
if (args.pix != None): pixels = args.pix
bunfile = args.file
bunOK = os.path.isfile(bunfile)
if not bunOK:
  print("file not found")
  sys.exit()

# Open the file for reading
print("reading ",bunfile)
data = sdds.SDDS(0)
data.load(bunfile)
if args.listpar:
    print()
    print("Parameters")
    print("==========")
    data.listParameters()
if args.listcol:
    print()
    print("Columns")
    print("=======")
    data.listColumns()

charge = 1.0e12 * data.getParameterValue("Charge")

x = np.array(data.getColumnData("x"))
DistProp(x, prt="x [m] : ")
xp = np.array(data.getColumnData("xp"))
DistProp(xp, prt="xp [] : ")
y = np.array(data.getColumnData("y"))
DistProp(y, prt="y [m] : ")
yp = np.array(data.getColumnData("yp"))
DistProp(yp, prt="yp [] : ")
p = np.array(data.getColumnData("p"))
t = np.array(data.getColumnData("t"))

(t0, tau) = DistProp(t, prt="arrival time t [s] : ")
(betagamma, prms) = DistProp(p, prt="particle momentum p [m_e*c] : ")
dt = t-t0

Np = len(x)
E = 0.511*np.sqrt(np.square(p)+1.0)
(E0, erms) = DistProp(E, prt="particle total energy E_tot [MeV] : ")
dE = E-E0

print("E-t correlation :")
MatXT = np.matrix([np.ones_like(1.0e12*dt),1.0e12*dt,np.square(1.0e12*dt)])
MatX = MatXT.transpose(1,0)
MatXTXI = np.dot(MatXT,MatX).getI()
MatXTY = np.dot(MatXT,E)[0]
corr = MatXTXI*MatXTY.transpose(1,0)
# print(corr)
print("%g keV/ps" % (1e3*corr[1,0]))
print("%g keV/ps²" % (1e3*corr[2,0]))

ex_rms = 1.0e6 * betagamma * np.sqrt((np.dot(x,x)*np.dot(xp,xp)-pow(np.dot(x,xp),2))/pow(float(Np),2))
ey_rms = 1.0e6 * betagamma * np.sqrt((np.dot(y,y)*np.dot(yp,yp)-pow(np.dot(y,yp),2))/pow(float(Np),2))
ez_rms = 1.0e15 * np.sqrt((np.dot(dE,dE)*np.dot(dt,dt)-pow(np.dot(dE,dt),2))/pow(float(Np),2))

string = r'p*c = %5.3f MeV' % E0 + r'  $\beta \gamma$ = %5.3f' % betagamma + '\n' + \
  '%d particles' % Np + '   Q=%.1f pC' % charge + '\n' + \
  r'$\sigma_E$ = %6.3f keV' % (1e3*erms) + '\n' + \
  r'$\sigma_t$ = %1.3f ps' % (1e12*tau) + '\n' + \
  r'$\epsilon^n_{x, RMS}$ = %2.3f $\mu m$' % ex_rms  + '\n' + \
  r'$\epsilon^n_{y, RMS}$ = %2.3f $\mu m$' % ey_rms  + '\n' + \
  r'$\epsilon_{z, RMS}$ = %3.1f keV ps' % ez_rms
fig = plt.figure(figsize=(16,10),dpi=80)
fig.text(0.68,0.93, string, color='b', size=18.0 , ha='left', va='top', linespacing=2, family='serif')
PlotPS(1000.0*x,1000.0*xp,"$x$ [mm]", "$p_x$ [mrad]", rect_dens = [0.03, 0.55, 0.3, 0.4])
PlotPS(1000.0*y,1000.0*yp,"$y$ [mm]", "$p_y$ [mrad]", rect_dens = [0.35, 0.55, 0.3, 0.4])
PlotPS(1000.0*x,1000.0*y,"$x$ [mm]", "$y$ [mm]", rect_dens = [0.03, 0.08, 0.3, 0.4], square=True)
PlotPS(1e12*dt,E,"$dt$ [ps]", "$E [MeV]$", rect_dens = [0.35, 0.08, 0.3, 0.4], center=False)

if args.currpro:
    hist, bins = np.histogram(1e12*dt, bins=100)
    delta = bins[1]-bins[0]
    center = (bins[:-1] + bins[1:]) / 2
    current = hist * np.fabs(charge)/delta/Np
    axc = plt.axes([0.68, 0.08, 0.3, 0.4])
    axc.bar(center, current, align='center', width=delta, zorder=10)
    axc.grid(True, zorder=5)
    axc.set_xlabel("$dt$ [ps]")
    axc.set_ylabel("$I$ [A]")
else:
    PlotPS(1e12*dt,1000*(E-c0-c1*1e12*dt-c2*np.square(1e12*dt)),"$dt$ [ps]", "$dE [keV]$", rect_dens = [0.68, 0.08, 0.3, 0.4], center=False)

if (args.img != None):
    plt.savefig(args.img)
else:
    plt.show()
