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
    \brief Radiation emitted by a point charge in a waveguided undulator

    @author Ulf Lehnert
    @date 27.8.2019
    @file waveguide_wiggler.cpp
    
    This test case tracks a single particle in an undulator field.
    The particle energy is 24 MeV it's charge 1 nC.
    
    That particle propagates through an undulator of 40 periods with 50 mm
    preiod length. The particle starts at z=0, the undulator is centered at z=2.0m.
    At z=4m the produced radiation is observed. A parallel plate waveguide is simulated
    by mirror particles of the tracked particle. With a sufficient number
    of mirror trajectories a waveguide up to the observation screen is obtained.
    
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
#include "screen.h"
#include "undulator.h"
#include <iostream>
#include <fstream>
#include <time.h>

int NOTS = 6000;	// number of time steps
int NOM = 40;       // number of mirror reflections

int main()
{
    double gamma = 14.0/0.511;
    double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    double betagamma = sqrt(gamma * gamma - 1.0);
    printf("beta =  %12.9g\n", beta);
    printf("gamma =  %12.9g\n", gamma);
    printf("c*p =  %12.9g MeV\n", 1e-6 * mecsquared * betagamma);

    double lambdar = 0.01/85.0;
    printf("Radiation Wavelength =  %6.3f µm\n", lambdar * 1.0e6);
    printf("Radiation Wavenumber =  %6.3f 1/cm\n", 0.01/lambdar);
    printf("Radiation Frequency =  %6.3g THz\n", SpeedOfLight/lambdar*1.0e-12);

    double lambda = 0.110;
    double N = 40;
    printf("Undulator Period = %9.6g m\n", lambda);
    printf("N = %9.6g\n", (double)N);
    double Ksq = lambdar/lambda * 2*gamma*gamma;
    double K = sqrt(2.0*(Ksq-1.0));
    double B = K / 93.4 / lambda;
    printf("B =  %9.6g T\n", B);
    PlanarUndulator* Undu = new PlanarUndulator(Vector(0.0, 0.0, 3.0));
    Undu->Setup(B, lambda, N);

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
    
    // Tracking should be done for 5.5 m in lab space corresponding to tau [s].
    // Inside the undulator we have an additional pathlength of one radiation
    // wavelength per period. The radiation wavelength already includes the
    // velocity of the particles. Outside the undulator the electron moves with beta*SpeedOfLight
    double tau = (double)N * (lambda + lambdar) / SpeedOfLight + (5.5 - (double)N * lambda) / (beta * SpeedOfLight);
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
	
	// compute the radiation observed on a finite screen 3.0m downstream the undulator center
	// this observer only sees the direct particle
	double t0 = (6.0 - 10*lambdar) / SpeedOfLight;
	double dt = lambdar / 12.0 / SpeedOfLight;
    ScreenObserver screenObs = ScreenObserver(
        "Waveguide_Screen_ObsRadField_direct.h5",
    	Vector(0.0, 0.0, 6.0),      // position
    	Vector(0.002, 0.0, 0.0),    // dx
    	Vector(0.0, 0.0002, 0.0),   // dy
    	101,                        // unsigned int nx,
    	51,                         // unsigned int ny,
    	t0,
    	dt,
    	6000);                      // NOTS
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
    
    // create a bunch including mirror particles
    ChargedParticle *e0 = new ChargedParticle(electron);
    Bunch *mirrored = new Bunch();
    mirrored->Add(e0);
	// the waveguide of 10mm height extends up to infinity
	// mirror planes are at y= +/- 5mm, 10mm, 15mm, ... with reversed polarity
	// even numbers of mirrors constitute a shift without polarity change
	for (int m=1; m<=NOM; m++)
	{
        ChargedParticle *e1 = new ChargedParticle(electron);
        ChargedParticle *e2 = new ChargedParticle(electron);
        if (m%2 == 0)
        {
            e1->Shift(Vector(0.0,0.010*m,0.0), 1.0);
            e2->Shift(Vector(0.0,-0.010*m,0.0), 1.0);
        }
        else
        {
            e1->Mirror(Vector(0.0,0.005*m,0.0), Vector(0.0,1.0,0.0), -1.0);
            e2->Mirror(Vector(0.0,-0.005*m,0.0), Vector(0.0,-1.0,0.0), -1.0);
        };
        mirrored->Add(e1);
        mirrored->Add(e2);
    }
	// compute the radiation observed on a finite screen 3.0m downstream the undulator center
	// this observer sees the mirror particles in addition
    ScreenObserver mirrorObs = ScreenObserver(
        "Waveguide_Screen_ObsRadField_with_mirrors.h5",
    	Vector(0.0, 0.0, 6.0),      // position
    	Vector(0.002, 0.0, 0.0),    // dx
    	Vector(0.0, 0.0002, 0.0),   // dy
    	101,                        // unsigned int nx,
    	51,                         // unsigned int ny,
    	t0,
    	dt,
    	4000);                      // NOTS
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);
    mirrorObs.integrate(mirrored);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop_time);
    printf("finished radiation calculation on screen.\n");
    elapsed = stop_time.tv_sec-start_time.tv_sec + 1e-9*(stop_time.tv_nsec-start_time.tv_nsec);
    printf("time elapsed : %9.3f s\n",elapsed);
    // write screen time traces
    try
    { 
    	mirrorObs.WriteTimeDomainFieldHDF5();
    	printf("Screen observer time domain field written - \033[1;32m OK\033[0m\n");
    }
    catch (exception& e) { cout << e.what() << endl;}
    
    // clean up
    delete lattice;
    // deleting the beam automatically
    // deletes all bunches and particles belonging to it
    delete beam;
    delete mirrored;
   
    return 0;
}
