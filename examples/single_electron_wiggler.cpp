/*=========================================================================

  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers

  Module:    undulator example (single particle)

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
    \brief Radiation emitted by a point charge in an undulator

    @author Ulf Lehnert
    @date 7.5.2018
    @file single_electron_wiggler.cpp
    
    This test case tracks a single particle in an undulator field.
    The particle energy is 24 MeV it's charge 1 nC.
    
    That particle propagates through an undulator of 8 periods with 300 mm
    preiod length. The particle starts at z=0, the undulator is centered at z=2.0m.
    At z=10m the produced radiation is observed.
    
    At c*t=2.0m a snapshot of the particle distribution and the local fields is generated.

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>

#include "bunch.h"
#include "beam.h"
#include "global.h"
#include "logger.h"
#include "observer.h"
#include "particle.h"
#include "fields.h"
#include "undulator.h"
#include <iostream>
#include <fstream>
#include <time.h>

int NOTS = 2000;	// number of time steps

int main()
{
    double B = 0.10315;
    double lambda = 0.300;
    double N = 8;
    PlanarUndulator* Undu = new PlanarUndulator(Vector(0.0, 0.0, 2.0));
    Undu->Setup(B, lambda, N);

    printf("B =  %9.6g T\n", B);
    printf("Undulator Period = %9.6g m\n", lambda);
    printf("N = %9.6g\n", (double)N);
    double gamma = 50.881;
    double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    double betagamma = sqrt(gamma * gamma - 1.0);
    double K = 0.934 * 100 * lambda * B;
    double lambdar = (lambda / (2 * gamma * gamma)) * (1 + K * K / 2);
    printf("beta =  %12.9g\n", beta);
    printf("gamma =  %12.9g\n", gamma);
    printf("c*p =  %12.9g MeV\n", 1e-6 * mecsquared * betagamma);
    printf("Radiation Wavelength =  %6.3f mm\n", lambdar * 1.0e3);
    printf("Radiation Frequency =  %6.3g THz\n", SpeedOfLight/lambdar*1.0e-12);
    
    // a simple lattice with just the Undulator Field
    Lattice* lattice = new Lattice;
    lattice->addElement(Undu);

    // setup a bunch with 
    // one single particle of 1nC
    // positioned at the origin (by default)
    double ch = 1.0e-9 / ElementaryCharge;
    ChargedParticle *electron = new ChargedParticle(-ch,ch);
    electron->initTrajectory(0.0,
        Vector(0.0,0.0,0.0),
        Vector(0.0,0.0,betagamma),
        Vector(0.0,0.0,0.0));
    Bunch *single = new Bunch();
    single->Add(electron);
    
    // Tracking should be done for 4 m in lab space corresponding to tau [s].
    // Inside the undulator we have an additional pathlength of one radiation
    // wavelength per period. The radiation wavelength already includes the
    // velocity of the particles. Outside the undulator the electron moves with beta*SpeedOfLight
    double tau = (double)N * (lambda + lambdar) / SpeedOfLight + (4.0 - (double)N * lambda) / (beta * SpeedOfLight);
    double deltaT = tau / NOTS;

    // track a beam of one bunche
    // containing only the single (reference) electron
    Beam *beam = new Beam();
    beam->Add(single);
    
    // setup for the tracking procedure
    beam->setTimeStep(deltaT);
    beam->InitVay(lattice);

    // do the tracking of the beam
    printf("tracking particles ...\n");
    fflush(stdout);
    // record the start time
    timespec start_time, stop_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);
    for (int step=0; step<NOTS; step++)
    {
    	beam->StepVay(lattice);
    }
    // record the finish time
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop_time);
    printf("finished tracking particles.\n");
    double elapsed = stop_time.tv_sec-start_time.tv_sec + 1e-9*(stop_time.tv_nsec-start_time.tv_nsec);
    printf("time elapsed during tracking : %9.3f s\n",elapsed);
	
    // compute the radiation of the single electron on axis
    double z0 = 2.0 + 10.0;
    double t0 = z0/SpeedOfLight - 1.0e-12;
    PointObserver singleObs = PointObserver(
	    "SingleParticle_ObsRadField.sdds", Vector(0.0, 0.0, z0), t0, 0.05e-13, 3000);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);
    singleObs.integrate(single);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop_time);
    printf("finished on-axis radiation calculation.\n");
    elapsed = stop_time.tv_sec-start_time.tv_sec + 1e-9*(stop_time.tv_nsec-start_time.tv_nsec);
    printf("time elapsed : %9.3f s\n",elapsed);
    // write field time traces
    try
    {
    	singleObs.WriteTimeDomainFieldSDDS();
	    printf("SDDS time domain field written - \033[1;32m OK\033[0m\n");
    }
    catch (exception& e) { cout << e.what() << endl;}
	
	// compute the radiation observed on a finite screen
    ScreenObserver screenObs = ScreenObserver(
        "SingleParticle_Screen_ObsRadField.h5",
    	Vector(0.0, 0.0, z0),       // position
    	Vector(0.004, 0.0, 0.0),    // dx
    	Vector(0.0, 0.004, 0.0),    // dy
    	201,                        // unsigned int nx,
    	201,                        // unsigned int ny,
    	t0,
    	0.5e-13,                    // double dt,
    	1000);                      // NOTS
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);
    screenObs.integrate(single);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop_time);
    printf("finished radiation calculation on screen.\n");
    elapsed = stop_time.tv_sec-start_time.tv_sec + 1e-9*(stop_time.tv_nsec-start_time.tv_nsec);
    printf("time elapsed : %9.3f s\n",elapsed);
    // write screen time traces
    try
    { 
    	screenObs.WriteTimeDomainFieldHDF5();
    	printf("Screen observer time domain field written - \033[1;32m OK\033[0m\n");
    }
    catch (exception& e) { cout << e.what() << endl;}

    // clean up
    delete lattice;
    // deleting the beam automatically
    // deletes all bunches and particles belonging to it
    delete beam;
    
    return 0;
}
