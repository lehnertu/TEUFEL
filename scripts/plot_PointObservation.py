#!/usr/bin/env python3
# coding=UTF-8

import sys, sdds, time
import os.path
import argparse
import numpy as np
import matplotlib.pyplot as plt

# magnetic field constant in N/A²
mu0 = 4*np.pi*1e-7

parser = argparse.ArgumentParser()
parser.add_argument('files', default=[], nargs='+', help='file name(s) of the SDDS file(s) with the observation field data')
parser.add_argument('--list_params', dest='listpar',
  action='store_const', const=True, default=False, help='list all parameters available in the file')
parser.add_argument('--list_columns', dest='listcol',
  action='store_const', const=True, default=False, help='list all columns available in the file')
parser.add_argument('--spectrum', dest='spectrum',
  action='store_const', const=True, default=False, help='show spectrum of horizontal polarization')
parser.add_argument('-fmax', help="plot range for the spectrum in Hz", dest="fmax", type=float)

args = parser.parse_args()

fig = plt.figure(1,figsize=(12,9))
left, width = 0.15, 0.80
rect1 = [left, 0.55, width, 0.40]  #left, bottom, width, height
rect2 = [left, 0.08, width, 0.40]
ax1 = fig.add_axes(rect1)
ax4 = fig.add_axes(rect2, sharex=ax1)

if args.spectrum:
    fig2 = plt.figure(2,figsize=(12,9))
    ax2 = fig2.add_axes([0.1, 0.1, 0.8, 0.8])

for file in args.files:

    fileOK = os.path.isfile(file)
    if not fileOK:
        print("file not found")
        sys.exit()

    print("reading ",file)
    data = sdds.SDDS(0)
    data.load(file)
    if args.listpar:
        data.listParameters()
    if args.listcol:
        data.listColumns()

    if ("t" in data.columnName):
        t = np.array(data.getColumnData("t"))*1e9
    if ("t0" in data.parameterName and
        "dt" in data.parameterName and
        "NumberTimeSteps" in data.parameterName):
        t0 = data.getParameterValue("t0")
        dt = data.getParameterValue("dt")
        nots = data.getParameterValue("NumberTimeSteps")
        print(r'%d steps starting at t0=%12.3g s step dt=%12.3g s' % (nots,t0,dt))
        t = (np.linspace(t0, t0+nots*dt, num=nots)+0.5*dt)*1e9
    else:
        print("unable to read the timing")
        sys.exit()

    Ex = np.array(data.getColumnData("Ex"))
    Ey = np.array(data.getColumnData("Ey"))
    Ez = np.array(data.getColumnData("Ez"))
    Bx = np.array(data.getColumnData("Bx"))
    By = np.array(data.getColumnData("By"))
    Bz = np.array(data.getColumnData("Bz"))
    
    EVec = np.array([Ex, Ey, Ez]).transpose()
    BVec = np.array([Bx, By, Bz]).transpose()
    # Poynting vector in V/m * (N/(A m)) / (N/A²) = W/m²
    SVec = np.cross(EVec, BVec) / mu0
    if 'dt' in globals():
        print('energy flow density = ', SVec.sum(axis=0)*dt, " Ws/m²")

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

    # compute and show spectrum for horizontal polarization
    if args.spectrum:
        nf = nots
        fmax = 1.0/dt
        f = np.linspace(0.0,fmax,nf)[:nf//2]
        df=1.0/dt/nf
        spectE = np.fft.fft(Ex)[:nf//2]
        spectB = np.fft.fft(By)[:nf//2]
        amplit = np.real(spectE*np.conj(spectB)/mu0)*2*dt/(df*nf)

        l21 = ax2.plot(1e-12*f, 1e6*1e12*amplit, "r-")
        if args.fmax != None:
            ax2.set_xlim(0.0,1e-12*args.fmax)
        ax2.set_xlabel(r'$f$ [THz]')
        ax2.set_ylabel(r'spectral power density on axis   $dI/dA/df$ [$\mu$J/(m$^2$ THz)]')
        ax2.yaxis.label.set_color('r')
        ax2.tick_params(axis='y', colors='r')

plt.show()
