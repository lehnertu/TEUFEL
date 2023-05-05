#!/usr/bin/env python3
# coding=UTF-8

import scipy.constants as sc
import os
import numpy as np
import h5py
import sdds

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
    def read(cls, filename, printout=True):
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
            if printout: print("Nx=%d Ny=%d Nots=%d" % (screen.Nx,screen.Ny,screen.nots))
            sd = np.array(sn)
            screen.origin = sd[0]
            if printout: print("origin= (%g, %g, %g) m" % (screen.origin[0],screen.origin[1],screen.origin[2]))
            screen.normal = sd[1]
            if printout: print("normal= (%g, %g, %g) m" % (screen.normal[0],screen.normal[1],screen.normal[2]))
            screen.dX = sd[2]
            screen.dx = np.linalg.norm(screen.dX)
            screen.ex = screen.dX / screen.dx
            if printout: print("dX= (%g, %g, %g) m" % (screen.dX[0],screen.dX[1],screen.dX[2]))
            screen.dY = sd[3]
            screen.dy = np.linalg.norm(screen.dY)
            screen.ey = screen.dY / screen.dy
            if printout: print("dY= (%g, %g, %g) m" % (screen.dY[0],screen.dY[1],screen.dY[2]))
            screen.t0 = sn.attrs.get('t0')
            screen.dt = sn.attrs.get('dt')
            screen.dtx = sn.attrs.get('dtx')
            screen.dty = sn.attrs.get('dty')
            if printout: print("t0=%g dt=%g dtx=%g dty=%g" % (screen.t0, screen.dt, screen.dtx, screen.dty))
        else:
            print("WARNING : old file format")
            pos = hdf['ObservationPosition']
            posdata = np.array(pos)
            screen.Nx = pos.attrs.get('Nx')
            screen.Ny = pos.attrs.get('Ny')
            screen.xcenter = screen.Nx//2
            screen.ycenter = screen.Ny//2
            screen.nots = field.attrs.get('NOTS')
            if printout: print("Nx=%d Ny=%d Nots=%d" % (screen.Nx,screen.Ny,screen.nots))
            screen.origin = posdata[screen.xcenter,screen.ycenter]
            if printout: print("origin= (%g, %g, %g) m" % (screen.origin[0],screen.origin[1],screen.origin[2]))
            screen.dX = posdata[screen.xcenter+1,screen.ycenter] - screen.origin
            screen.dx = np.linalg.norm(screen.dX)
            screen.ex = screen.dX / screen.dx
            if printout: print("dX= (%g, %g, %g) m" % (screen.dX[0],screen.dX[1],screen.dX[2]))
            screen.dY = posdata[screen.xcenter,screen.ycenter+1] - screen.origin
            screen.dy = np.linalg.norm(screen.dY)
            screen.ey = screen.dY / screen.dy
            if printout: print("dY= (%g, %g, %g) m" % (screen.dY[0],screen.dY[1],screen.dY[2]))
            screen.normal = np.cross(screen.dX,screen.dY)*-1.0
            screen.normal = screen.normal / np.linalg.norm(screen.normal)
            if printout: print("normal= (%g, %g, %g) m" % (screen.normal[0],screen.normal[1],screen.normal[2]))
            screen.t0 = field.attrs.get('t0')
            screen.dt = field.attrs.get('dt')
            screen.dtx = 0.0
            screen.dty = 0.0
            if printout: print("t0=%g dt=%g dtx=%g dty=%g" % (screen.t0, screen.dt, screen.dtx, screen.dty))
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

    def df(self):
        """
        Report the frequency resolution of the recorded field traces
        """
        return 1.0/self.dt/self.nots
        
    def frequency(self,index):
        """
        get the frequency of the queried spectral channel
        """
        return index/self.dt/self.nots
        
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
        SVec = np.cross(EVec, BVec) / sc.mu_0
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
        amplit = np.abs(spectE*np.conj(spectB)/sc.mu_0)*2*self.dt/(df*nf)
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
        amplit = np.abs(spectE*np.conj(spectB)/sc.mu_0)*2*self.dt/(df*nf)
        return (f, amplit)

# =============================================================================

class TeufelWatch():
    """
    This class describes a watch point of particle coordinates
    written by TEUFEL during particle tracking.
    """
    
    def __init__(self):
        """
        Create an empty TeufelWatch object.
        """
        self.NOP = 0

    @classmethod
    def read(cls, filename):
        """
        Import a TeufelWatch object from a file.
        Both SDDS and HDF5 file formats are handled.
        """
        watch = cls.__new__(cls)
        # determine the file format
        f = open(filename, "rb")
        bytes = f.read(4)
        is_sdds = (bytes == b'SDDS')
        is_hdf5 = (bytes == b'\x89HDF')
        f.close()
        if not (is_sdds or is_hdf5):
            raise ValueError("File is neither SDDS nor HDF5 format")
        # read HDF5 data
        if is_hdf5:
            print("reading ",filename)
            hdf = h5py.File(filename, "r")
            electrons = hdf['electrons']
            a = np.array(electrons)
            hdf.close()
            data = a.transpose()
            watch.x = data[0]
            watch.y = data[1]
            watch.z = data[2]
            watch.bgx = data[3]
            watch.bgy = data[4]
            watch.bgz = data[5]
        # read SDDS data
        if is_sdds:
            print("reading ",filename)
            data = sdds.SDDS(0)
            data.load(filename)
            watch.x = np.array(data.getColumnData("x"))
            xp = np.array(data.getColumnData("xp"))
            watch.y = np.array(data.getColumnData("y"))
            yp = np.array(data.getColumnData("yp"))
            p = np.array(data.getColumnData("p"))
            watch.z = np.array(data.getColumnData("z"))
            watch.bgz = p*np.sqrt(1.0-np.square(xp)-np.square(yp))
            watch.bgx = bgz*xp
            watch.bgy = bgz*yp
        watch.NOP = watch.x.shape[0]
        print(f'have read {watch.NOP} particles.')
        # compute more particle properties
        watch.gamma = np.sqrt(np.square(watch.bgx)+np.square(watch.bgy)+np.square(watch.bgz)+1.0)
        # object is ready
        return watch

    def parameters(self):
        """
        This method computes a number of parameters of the beam
        and returns a dictionary of these values.
        N         number of particles
        gamma     average relativistic factor
        E_kin     average kinetic energy [eV]
        E_rms     energy spread [eV]
        x_rms     rms beam size x [m]
        y_rms     rms beam size y [m]
        z_rms     rms beam size z [m]
        t_rms     rms bunch length [s]
        ex_n_rms  normalized transverse emittance [m rad]
        ey_n_rms  normalized transverse emittance [m rad]
        ez_rms    longitudinal emittance [eV s]
        """
        mec2 = sc.m_e*sc.c**2/sc.e
        # per-particle computations
        bg = np.sqrt(np.square(self.bgx) + np.square(self.bgy) + np.square(self.bgz))
        gamma = np.sqrt(np.square(self.bgx) + np.square(self.bgy) + np.square(self.bgz) + 1.0)
        E_kin = (gamma-1.0)*mec2
        xp=np.arctan2(self.bgx,self.bgz)
        yp=np.arctan2(self.bgy,self.bgz)
        # mean values
        xmean = self.x.mean()
        ymean = self.y.mean()
        xpmean = xp.mean()
        ypmean = yp.mean()
        zmean = self.z.mean()
        bgmean = bg.mean()
        meanE = np.mean(E_kin)
        # rms values
        dx = self.x - xmean
        dy = self.y - ymean
        dxp = xp - xpmean
        dyp = yp - ypmean
        dz = self.z - zmean
        dt = dz/sc.c
        dE = E_kin - meanE
        xrms = np.sqrt(np.dot(dx,dx)/self.NOP)
        yrms = np.sqrt(np.dot(dy,dy)/self.NOP)
        zrms = np.sqrt(np.dot(dz,dz)/self.NOP)
        trms = np.sqrt(np.dot(dt,dt)/self.NOP)
        Erms = np.sqrt(np.dot(dE,dE)/self.NOP)
        ex = bgmean * np.sqrt((np.dot(dx,dx)*np.dot(dxp,dxp)-pow(np.dot(dx,dxp),2)))/self.NOP
        ey = bgmean * np.sqrt((np.dot(dy,dy)*np.dot(dyp,dyp)-pow(np.dot(dy,dyp),2)))/self.NOP
        ez = np.sqrt((np.dot(dE,dE)*np.dot(dt,dt)-pow(np.dot(dE,dt),2)))/self.NOP
        return {'N': self.NOP, 'x_mean':xmean, 'y_mean':ymean, 'z_mean':zmean, 'xp_mean':xp.mean(), 'yp_mean':yp.mean(),
            'bgx_mean':self.bgx.mean(), 'bgy_mean':self.bgy.mean(), 'bgz_mean':self.bgz.mean(),
            'E_kin':meanE, 'E_tot':gamma.mean()*mec2, 'gamma':gamma.mean(), 'betagamma':bgmean,
            'x_rms':xrms, 'y_rms':yrms, 'z_rms':zrms, 't_rms':trms, 'E_rms':Erms,
            'ex_n_rms':ex, 'ey_n_rms':ey, 'ez_rms':ez,}

# =============================================================================

def readTeufelBeamLog(filename):
    """
    This procedure reads a log of beam parameters (SDDS file)
    written by TEUFEL during particle tracking.
    It returns a dictionary of properties each containing a 
    list of values observed along the trajectory of the beam.
    """
    if not os.path.exists(filename):
        print('file not found.')
        return {}
    data = sdds.SDDS(0)
    data.load(filename)
    x = np.array(data.getColumnData("x_av"))
    y = np.array(data.getColumnData("y_av"))
    z = np.array(data.getColumnData("z_av"))
    if ("t" in data.columnName):
        t = np.array(data.getColumnData("t"))
    else:
        t = np.zeros_like(x);
    bgx = np.array(data.getColumnData("bgx_av"))
    bgy = np.array(data.getColumnData("bgy_av"))
    bgz = np.array(data.getColumnData("bgz_av"))
    if ("gamma" in data.columnName):
        gamma = np.array(data.getColumnData("gamma"))
    else:
        gamma = np.zeros_like(x);
    if ("delta" in data.columnName):
        delta = np.array(data.getColumnData("delta"))
    else:
        delta = np.zeros_like(x);
    if ("BF" in data.columnName):
        bf = np.array(data.getColumnData("BF"))
    else:
        bf = np.zeros_like(x);
    x_rms = np.array(data.getColumnData("x_rms"))
    y_rms = np.array(data.getColumnData("y_rms"))
    z_rms = np.array(data.getColumnData("z_rms"))
    return {'x_avg':x, 'y_avg':y, 'z_avg':z, 't':t, 'x_rms':x_rms, 'y_rms':y_rms, 'z_rms':z_rms,
        'bgx_avg':bgx, 'bgy_avg':bgy, 'bgz_avg':bgz, 'gamma':gamma, 'delta':delta, 'BF':bf}

