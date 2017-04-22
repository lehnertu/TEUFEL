/*=========================================================================

  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers

  Module:    undulator test case

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
    \brief Radiation from the ELBE U300 THz source

    @author Ulf Lehnert
    @date 12.4.2017
    @file elbe-u300.cpp
    
    This test case tracks a single particle in an undulator field.
    The particle energy is 24 MeV it corresponds to 4.37e8 electrons
    that is a charge of 70 pC.
    
    That particle propagates through an undulator of 8 periods with 300 mm
    preiod length. The particle starts at z=0, the undulator is centered at z=2.0m.
    At z=10m the produced radiation is observed.
    
    The program generates a trajectory dump elbe-u300_trajectory.sdds which
    can be used to plot the elctron trajectory.
    
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "observer.h"
#include "particle.h"
#include "fields.h"
#include "undulator.h"
#include <iostream>
#include <fstream>

int NOTS = 4000;    // number of time steps

int main()
{
    double B = 0.100;
    double lambda = 0.300;
    double N = 8;
    PlanarUndulator* Undu = new PlanarUndulator(Vector(0.0, 0.0, 2.0));
    Undu->Setup(B, lambda, N);

    printf("B =  %9.6g T\n", B);
    printf("Undulator Period = %9.6g m\n ", lambda);
    printf("N = %9.6g\n ", (double)N);
    double gamma = 48;
    double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    double betagamma = sqrt(gamma * gamma - 1.0);
    double K = 0.934 * 100 * lambda * B;
    double lambdar = (lambda / (2 * gamma * gamma)) * (1 + K * K / 2);
    printf("beta =  %12.9g\n", beta);
    printf("gamma =  %12.9g\n", gamma);
    printf("c*p =  %12.9g MeV\n", 1e-6 * mecsquared * betagamma);
    printf("Radiation Wavelength =  %6.3f mm\n", lambdar * 1.0e3);

    // a simple lattice with just the Undulator Field
    Lattice* lattice = new Lattice;
    lattice->addElement(Undu);

    // one single particle corresponding to 4.37e8 electrons
    ChargedParticle* electron = new ChargedParticle(-4.37e8,4.37e8);

    // initial position at the origin
    Vector X0 = Vector(0.0, 0.0, 0.0);
    // initial momentum of the particle
    Vector P0 = Vector(0.0, 0.0, betagamma);

    // Track the particle for 3.0 m in lab space.
    // Inside the undulator we have an additional pathlength of one radiation
    // wavelength per period. The radiation wavelength already includes the
    // velocity of the particles. Outside the undulator the electron moves with
    // beta*SpeedOfLight
    double tau = (double)N * (lambda + lambdar) / SpeedOfLight + (4.0 - (double)N * lambda) / (beta * SpeedOfLight);
    double deltaT = tau / NOTS;
    electron->TrackVay(NOTS, deltaT, X0, P0, lattice);

    // create a trajectory dump
    if (0 != electron->WriteSDDS("elbe-u300_Trajectory.sdds"))
    {
        printf("SDDS write \033[1;31m failed!\033[0m\n");
    }
    else
    {
        printf("SDDS file written - \033[1;32m OK\033[0m\n");
    }

    // compute the radiation on axis
    PointObserver Obs = PointObserver(Vector(0.0, 0.0, 10.0));
    Obs.GetTimeDomainTrace(electron);
    // dump it to a file
    if (0 != Obs.WriteTimeTraceSDDS("elbe-u300_RadTrace.sdds"))
    {
	printf("SDDS write \033[1;31m failed!\033[0m\n");
    }
    else
    {
	printf("SDDS file written - \033[1;32m OK\033[0m\n");
    }

    // clean up
    delete lattice;
    delete electron;

    return 0;
}
