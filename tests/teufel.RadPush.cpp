/*=========================================================================

  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers

  Module:    Homogeneous magnetic field test case

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
    \brief Radiation push of a free electron

    @author Ulf Lehnert
    @date 10.4.2019
    @file teufel.RadPush.cpp
    
    This test case tracks a single electron in a gaussian electromagnetic wave.
    The wavelength is 1 micron, the field strength choosen such that \f$ a_0 = 1 \f$.
    @li The particle should oscillate with an amplitude of \f$ a_0 / k \f$.
    @li The particle should be accelerated to \f$ \beta \gamma = 1 \f$ in the field direction.
    @li Due to the magnetic field component it should move in positive z direction
    with a peak momentum of \f$ \beta \gamma = a_0^2 / 2 \f$
    <br>
    
    For plotting the trajectory a log file rad_push_log.sdds is created.
    <br>
    
    @return The number of errors encountered in the above list of checks is reported.
    
 */

#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "logger.h"
#include "particle.h"
#include "fields.h"
#include "wave.h"

int NOTS = 2000;                // number of time steps

int main ()
{

    printf("\nTEUFEL - radiation push testcase\n");

	GaussianWave *wave = new GaussianWave();
	wave->Setup(1.0e-6, complex<double>(3.2107e12,0.0), 1.0);
	ElMagField F = wave->Field(0.0, Vector(0.0, 0.0, 0.0));
    double Ex = F.E().x;
    double I0 = (Ex*Ex)/2.0*EpsNull*SpeedOfLight;
    printf("I0 =  %9.6g W/cm²\n",1e-4*I0);
    double lambda = wave->getWavelength();
    printf("lambda =  %9.6g µm\n",1e6*lambda);
    double omega = 2*Pi*SpeedOfLight / wave->getWavelength();
    printf("omega =  %9.6g Hz\n",omega);
	double a0 = Ex / (mecsquared/SpeedOfLight*omega);
    double a0_k = SpeedOfLight*a0/omega;
    printf("a0 = %9.6g    a0/k = %9.6g\n",a0,a0_k);
    
    // a simple lattice with just the inpinging wave
    Lattice *lattice = new Lattice;
    lattice->addElement(wave);
    
    // one single electron
    ChargedParticle *electron = new ChargedParticle();
    Bunch *bunch = new Bunch();
    bunch->Add(electron);
    
    // initial position at the origin
    Vector X0 = Vector(0.0, 0.0, 0.0);
    // initial momentum of the particle
    Vector P0 = Vector(0.0, 0.0, 0.0);
    Vector A0 = Vector(0.0, 0.0, 0.0);
    electron->initTrajectory(0.0, X0, P0, A0);
    
    // count the errors
    int errors = 0;
    
    // Track the particle for 10 cycles of the wave
    double tau = wave->getWavelength() / SpeedOfLight;
    printf("simulation running %9.6g s\n",10.0*tau);
    double deltaT = 10*tau / NOTS;
    bunch->InitVay(deltaT, lattice);
    TrackingLogger<Bunch> *log = new TrackingLogger<Bunch>(bunch, "rad_push_log.sdds");
    log->update();
    
    // record maximum displacement in x-direction
    double xmin = 0.0;
    double xmax = 0.0;
    double bgxmin = 0.0;
    double bgxmax = 0.0;
    double bgzmax = 0.0;
    for (int i=0; i<NOTS; i++)
    {
    	bunch->StepVay(lattice);
    	log->update();
    	Vector x = electron->getPosition();
    	if (x.x<xmin) xmin = x.x;
    	if (x.x>xmax) xmax = x.x;
    	Vector bg = electron->getMomentum();
    	if (bg.x<bgxmin) bgxmin = bg.x;
    	if (bg.x>bgxmax) bgxmax = bg.x;
    	if (bg.z>bgzmax) bgzmax = bg.z;
    };
    double t = electron->getTime();
    Vector XP = electron->getPosition();
    printf("x(%9.6g s) =  (%9.6g,%9.6g,%9.6g) m\n",t,XP.x,XP.y,XP.z);
    double xamp = xmax-xmin;

    if (fabs(xamp-2*a0_k)/a0_k > 1e-3) {
		errors++;
		printf("error in oscillation amplitude %9.6g m - \033[1;31m test failed!\033[0m\n", xamp);
    } else {
		printf("oscillation amplitude %9.6g m - \033[1;32m OK\033[0m\n", xamp);
    }

    if ( fabs(bgxmin+1)>1e-3 || fabs(bgxmax-1)>1e-3) {
		errors++;
		printf("error in beta*gamma : (%9.6g, %9.6g) - \033[1;31m test failed!\033[0m\n", bgxmin, bgxmax);
    } else {
		printf("beta*gamma : (%9.6g, %9.6g) - \033[1;32m OK\033[0m\n", bgxmin,bgxmax);
    }

    if ( fabs(bgzmax-0.5)>1e-3 ) {
		errors++;
		printf("error in peak beta*gamma = %9.6g - \033[1;31m test failed!\033[0m\n", bgzmax);
    } else {
		printf("peak beta*gamma = %9.6g - \033[1;32m OK\033[0m\n", bgzmax);
    }

	int res = log->WriteBeamParametersSDDS();
	errors += res;
	
    // clean up
    delete lattice;
    delete electron;

    return errors;
}
