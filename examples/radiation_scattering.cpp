/*=========================================================================

  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers

  Module:    radiation scattering example

  Copyright (c) 2017 U. Lehnert

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/

/*!
    \brief Radiation scattering

    @author Ulf Lehnert
    @date 8.8.2019
    @file radiation_scattering.cpp
    
    This example creates a wavepacket (as a lattice element).
    The radiation field is observed on a screen 1m downstream.

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>

#include "fields.h"
#include "global.h"
#include "wave.h"
#include "observer.h"
#include <iostream>
#include <fstream>
#include <time.h>

int NOTS = 2000;	// number of time steps

int main()
{
    // this must be globally defined even though it is only used in the MPI version
    teufel::rank = 0;

    // wave packet with waist of the gaussian envelope at the origin
    GaussianWavePacket* src = new GaussianWavePacket(Vector(0.0,0.0,0.0));
    // wavelength 100 micron
    // amplitude 3.5 MV/m
    // Rayleigh range 1m
    // pulse duration 1 ps
    // start time 0 at the waist (origin)
    src->Setup(1.0e-4, complex<double>(3.5e6, 0.0), 1.0, 1.0e-12, 0.0);

    // a simple lattice with just the wave packet
    Lattice* lattice = new Lattice;
    lattice->addElement(src);

    // define one observation screen at the origin
    ScreenObserver screenObs = ScreenObserver(
        "scattering_Screen_origin.h5",
    	Vector(0.0, 0.0, 0.0),		// position
    	Vector(0.001, 0.0, 0.0),		// dx
    	Vector(0.0, 0.001, 0.0),		// dy
    	41,				// unsigned int nx,
    	41,				// unsigned int ny,
    	-5.0e-12,
    	1.0e-14,			// double dt,
    	1000);				// NOTS
    screenObs.integrate(lattice);
    try
    { 
    	screenObs.WriteTimeDomainFieldHDF5();
    	printf("Screen observer time domain field written - \033[1;32m OK\033[0m\n");
    }
    catch (exception& e) { cout << e.what() << endl;}

    // define one observation screen at the origin
    double z0 = 1.0;
    double t0 = z0/SpeedOfLight - 5.0e-12;
    screenObs = ScreenObserver(
        "scattering_Screen_1m.h5",
    	Vector(0.0, 0.0, z0),		// position
    	Vector(0.001, 0.0, 0.0),		// dx
    	Vector(0.0, 0.001, 0.0),		// dy
    	41,				// unsigned int nx,
    	41,				// unsigned int ny,
    	t0,
    	1.0e-14,			// double dt,
    	1000);				// NOTS
    screenObs.integrate(lattice);
    try
    { 
    	screenObs.WriteTimeDomainFieldHDF5();
    	printf("Screen observer time domain field written - \033[1;32m OK\033[0m\n");
    }
    catch (exception& e) { cout << e.what() << endl;}

    // clean up (deleting the lattice automatically deletes all elements)
    delete lattice;

    return 0;
}
