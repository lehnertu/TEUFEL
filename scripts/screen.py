#!/usr/bin/env python3
# coding=UTF-8

import sys, time
import os.path
import numpy as np
import h5py

mu0 = 4*np.pi*1e-7

class TeufelScreen():
    """
    This class describes a time-domain optical field recorded
    on a screen during a TEUFEL run. It is used to hold the data
    for post-processing.
    """
    
    def __init__(screen):
        """
        Create an empty TeufelScreen object.
        """

    @classmethod
    def read(cls, filename):
        """
        Import a TeufelScreen object from a file.
        """
        screen = cls()
        hdf = h5py.File(filename, "r")
        # get the field data
        field = hdf['ElMagField']
        screen.A = np.array(field)
        # get the geomtry information
        sn = hdf['Screen']
        screen.Nx = sn.attrs.get('Nx')
        screen.Ny = sn.attrs.get('Ny')
        screen.nots = sn.attrs.get('NOTS')
        print("Nx=%d Ny=%d Nots=%d" % (screen.Nx,screen.Ny,screen.nots))
        sd = np.array(sn)
        screen.origin = sd[0]
        print("origin= (%g, %g, %g) m" % (screen.origin[0],screen.origin[1],screen.origin[2]))
        screen.normal = sd[1]
        print("normal= (%g, %g, %g) m" % (screen.normal[0],screen.normal[1],screen.normal[2]))
        screen.dX = sd[2]
        screen.dx = np.linalg.norm(screen.dX)
        screen.ex = screen.dX / screen.dx
        print("dX= (%g, %g, %g) m" % (screen.dX[0],screen.dX[1],screen.dX[2]))
        screen.dY = sd[3]
        screen.dy = np.linalg.norm(screen.dY)
        screen.ey = screen.dY / screen.dy
        print("dY= (%g, %g, %g) m" % (screen.dY[0],screen.dY[1],screen.dY[2]))
        screen.t0 = sn.attrs.get('t0')
        screen.dt = sn.attrs.get('dt')
        screen.dtx = sn.attrs.get('dtx')
        screen.dty = sn.attrs.get('dty')
        print("t0=%g dt=%g dtx=%g dty=%g" % (screen.t0, screen.dt, screen.dtx, screen.dty))
        hdf.close()
        # done with the file
        screen.xcenter = screen.Nx//2
        screen.ycenter = screen.Ny//2
        # object is ready
        return screen
        
    def screenPos(self,ix,iy):
        """
        get the position [m] of an indexed pixel on the screen
        relative to its origin
        """
        x = self.dx * (ix-self.xcenter)
        y = self.dy * (iy-self.ycenter)
        return (x,y)

    def getFieldTrace(self,ix,iy):
        """
        get the fields recorded in one pixel
        E [V/m] and B [T] separately as arrays of vectors.
        """
        trace = self.A[ix,iy]
        data = trace.transpose()
        Ex = data[0]
        Ey = data[1]
        Ez = data[2]
        Bx = data[3]
        By = data[4]
        Bz = data[5]
        EVec = np.array([Ex, Ey, Ez]).transpose()
        BVec = np.array([Bx, By, Bz]).transpose()
        return (EVec, BVec)
        
    def totalPowerDensity(self,ix,iy):
        """
        get the total incident radiation energy [J/m²] per unit area of the screen
       
        """
        (EVec, BVec) = self.getFieldTrace(ix,iy)
        SVec = np.cross(EVec, BVec) / mu0
        return np.dot ( SVec.sum(axis=0), -self.normal) * self.dt
        
    def getInPlaneFields(self,ix,iy):
        """
        get the in-plane field components as recorded in one pixel
        """
        (EVec, BVec) = self.getFieldTrace(ix,iy)
        EH = np.dot (EVec, self.ex)
        EV = np.dot (EVec, self.ey)
        BH = np.dot (BVec, self.ex)
        BV = np.dot (BVec, self.ey)
        return (EH, EV, BH, BV)

    def spectralDensityH(self,ix,iy):
        """
        get the spectral power density [J/(m² THz)] of the horizontal polarization
        """  
        nf = self.nots
        fmax = 1.0/self.dt
        f = np.linspace(0.0,fmax,nf)[:nf//2]
        df=1.0/self.dt/nf
        (EH, EV, BH, BV) = self.getInPlaneFields(ix,iy)
        spectE = np.fft.fft(EH)[:nf//2]
        spectB = np.fft.fft(BV)[:nf//2]
        amplit = np.real(spectE*np.conj(spectB)/mu0)*2*self.dt/(df*nf)
        return (f, amplit)

