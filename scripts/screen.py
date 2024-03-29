#!/usr/bin/env python3
# coding=UTF-8

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
        if hdf.get('Screen', getclass=True) == h5py._hl.dataset.Dataset:
            # new version of the file format
            sn = hdf['Screen']
            screen.Nx = sn.attrs.get('Nx')
            screen.Ny = sn.attrs.get('Ny')
            screen.xcenter = screen.Nx//2
            screen.ycenter = screen.Ny//2
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
        else:
            print("WARNING : old file format")
            pos = hdf['ObservationPosition']
            posdata = np.array(pos)
            screen.Nx = pos.attrs.get('Nx')
            screen.Ny = pos.attrs.get('Ny')
            screen.xcenter = screen.Nx//2
            screen.ycenter = screen.Ny//2
            screen.nots = field.attrs.get('NOTS')
            print("Nx=%d Ny=%d Nots=%d" % (screen.Nx,screen.Ny,screen.nots))
            screen.origin = posdata[screen.xcenter,screen.ycenter]
            print("origin= (%g, %g, %g) m" % (screen.origin[0],screen.origin[1],screen.origin[2]))
            screen.dX = posdata[screen.xcenter+1,screen.ycenter] - screen.origin
            screen.dx = np.linalg.norm(screen.dX)
            screen.ex = screen.dX / screen.dx
            print("dX= (%g, %g, %g) m" % (screen.dX[0],screen.dX[1],screen.dX[2]))
            screen.dY = posdata[screen.xcenter,screen.ycenter+1] - screen.origin
            screen.dy = np.linalg.norm(screen.dY)
            screen.ey = screen.dY / screen.dy
            print("dY= (%g, %g, %g) m" % (screen.dY[0],screen.dY[1],screen.dY[2]))
            screen.normal = np.cross(screen.dX,screen.dY)*-1.0
            screen.normal = screen.normal / np.linalg.norm(screen.normal)
            print("normal= (%g, %g, %g) m" % (screen.normal[0],screen.normal[1],screen.normal[2]))
            screen.t0 = field.attrs.get('t0')
            screen.dt = field.attrs.get('dt')
            screen.dtx = 0.0
            screen.dty = 0.0
            print("t0=%g dt=%g dtx=%g dty=%g" % (screen.t0, screen.dt, screen.dtx, screen.dty))
        # done with the file
        hdf.close()
        # object is ready
        return screen

    def write(self, filename):
        """
        Export a TeufelScreen object to an HDF5 file.
        """
        hf = h5py.File(filename, 'w')
        # data for the 'Screen' dataset
        sd = np.array([self.origin, self.normal, self.dX, self.dY])
        dataset = hf.create_dataset('Screen', data=sd)
        # add attributes for the geometry and sampling
        dataset.attrs['Nx'] = self.Nx
        dataset.attrs['Ny'] = self.Ny
        dataset.attrs['NOTS'] = self.nots
        dataset.attrs['t0'] = self.t0
        dataset.attrs['dt'] = self.dt
        dataset.attrs['dtx'] = self.dtx
        dataset.attrs['dty'] = self.dty
        # data for the 'ElMagField' dataset
        hf.create_dataset('ElMagField', data=self.A)
        # write to disk
        hf.close()
        
    def shape(self):
        """
        Report the shape of the screen.
        """
        return (self.Nx, self.Ny)

    def screenPos(self,ix,iy):
        """
        get the position [m] of an indexed pixel on the screen
        relative to its origin
        """
        x = self.dx * (ix-self.xcenter)
        y = self.dy * (iy-self.ycenter)
        return (x,y)

    def getStartTime(self,ix,iy):
        """
        get the start time [s] of the trace for an indexed pixel on the screen
        """
        return self.t0 + self.dtx*(ix-self.xcenter) + self.dty*(iy-self.ycenter)

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
    
    def setFieldTrace(self,ix,iy,EVec,BVec):
        """
        eet the fields for one pixel from 
        E [V/m] and B [T] arrays of vectors.
        """
        ET = EVec.transpose()
        Ex = ET[0]
        Ey = ET[1]
        Ez = ET[2]
        BT = BVec.transpose()
        Bx = BT[0]
        By = BT[1]
        Bz = BT[2]
        trace = np.array([Ex, Ey, Ez, Bx, By, Bz]).transpose()
        self.A[ix,iy] = trace
    
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
        amplit = np.abs(spectE*np.conj(spectB)/mu0)*2*self.dt/(df*nf)
        return (f, amplit)

    def spectralDensityV(self,ix,iy):
        """
        get the spectral power density [J/(m² THz)] of the vertical polarization
        """  
        nf = self.nots
        fmax = 1.0/self.dt
        f = np.linspace(0.0,fmax,nf)[:nf//2]
        df=1.0/self.dt/nf
        (EH, EV, BH, BV) = self.getInPlaneFields(ix,iy)
        spectE = np.fft.fft(EV)[:nf//2]
        spectB = np.fft.fft(BH)[:nf//2]
        amplit = np.abs(spectE*np.conj(spectB)/mu0)*2*self.dt/(df*nf)
        return (f, amplit)

