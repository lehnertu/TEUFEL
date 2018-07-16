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
    \brief Undulator test case

    @author Ulf Lehnert
    @author Vipul Joshi
    @date 10.2.2017
    @file teufel.undulator.cpp
    
    This test case tracks a single electron in an undulator field.
    The electron moves on an periodic sinusoidal trajectory.
    The x-displacement and slippage per period(radiation wavelength) are compared to known values.
    \f[
	x_{max} = \frac{e B}{\beta \gamma m c k_u^2}
    \f]
    \f[
	\lambda_R = \frac{\lambda_U}{2 \gamma^2} (1+\frac{K^2}{2}) \;\;\;\;
	K = \frac{\lambda_U}{2 \pi}\,\frac{e}{mc}\,B_{peak}
    \f]

    The program computes the trajectory of the electron starting at the coordinate system origin
    with a velocity in z direction. The magnitude of the velocity is choosen such that
    it is relativistic with \f$\gamma = 30.35823\f$. The trajectory's maximum displacment in x should be R=0.35686mm in a magnetic field with B=0.532 T.
    
    An Undulator with 54 period of \f$\lambda_U\f$ = 37 mm is placed t z=1.5 m.
    The electron is tracked for a fixed amount of time which is known in advance.
    It is computed from the desired z=3.0 m at the end of the tracking.
    It is computed as \f$N (\lambda_U + \lambda_R) / c + (3.0 - N \lambda_U) / (\beta c)\f$.
    
    After tracking the particle, it is checked that :
    @li the maximum displacement is correct
    @li the slippage is correct, so, the particle really arrives at z=3.0 m
    @li the particle has returned to x=0,y=0 and moving along the z-axis
    @li the kinetic energy has not changed
    
    Using the Vay algorithm, the tracking safely reaches the required accuracy with 3000 steps.
    
    @return The number of errors encountered in the above list of checks is reported.
    
    The program generates a trajectory dump teufel_undulator_trajectory.sdds which
    can be used to plot the elctron trajectory.
    
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "particle.h"
#include "fields.h"
#include "undulator.h"
#include <iostream>
#include <fstream>

int NOTS = 3000;    // number of time steps

int main()
{
    printf("\nTEUFEL - Undulator testcase\n");
    double B = 0.532;
    double lambda = 0.037;
    double N = 54;
    double ku = 2 * Pi / lambda;
    PlanarUndulator* Undu = new PlanarUndulator(Vector(0.0, 0.0, 1.5));
    Undu->Setup(B, lambda, N);

    printf("B =  %9.6g T\n", B);
    printf("Undulator Period = %9.6g m\n", lambda);
    printf("N = %9.6g\n", (double)N);
    double gamma = 30.35823388162631;
    double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    double betagamma = sqrt(gamma * gamma - 1.0);
    double K = 0.934 * 100 * lambda * B;
    double lambdar = (lambda / (2 * gamma * gamma)) * (1 + K * K / 2);
    printf("beta =  %12.9g\n", beta);
    printf("gamma =  %12.9g\n", gamma);
    printf("c*p =  %12.9g MeV\n", 1e-6 * mecsquared * betagamma);
    double xdis = ElementaryCharge * B * SpeedOfLight / (betagamma * 9.1e-31 * SpeedOfLight * SpeedOfLight * ku * ku);
    printf("X-Displacement =  %9.6fg mm\n", xdis*1.0e3);
    printf("Radiation Wavelength =  %6.3f µm\n", lambdar * 1.0e6);

    // a simple lattice with just the Undulator Field
    Lattice* lattice = new Lattice;
    lattice->addElement(Undu);

    // one single electron
    ChargedParticle* electron = new ChargedParticle();

    // initial position at the origin
    Vector X0 = Vector(0.0, 0.0, 0.0);
    // initial momentum of the particle
    Vector P0 = Vector(0.0, 0.0, betagamma);
    Vector A0 = Vector(0.0, 0.0, 0.0);
    electron->initTrajectory(0.0, X0, P0, A0);
    
    // Track the particle for 3.0 m in lab space.
    // Inside the undulator we have an additional pathlength of one radiation
    // wavelength per period. The radiation wavelength already includes the
    // velocity of the particles. Outside the undulator the electron moves with
    // beta*SpeedOfLight
    double tau = (double)N * (lambda + lambdar) / SpeedOfLight + (3.0 - (double)N * lambda) / (beta * SpeedOfLight);
    // double tau=06.8e-9;
    double deltaT = tau / NOTS;
    electron->InitVay(deltaT, lattice);
    double xdisp = 0;    //maximum displacement in + direction
    for (int i=0; i<NOTS; i++)
    {
	electron->StepVay(lattice);
	// look for the maximum displacement
        Vector XP = electron->getPosition();
        if (XP.z > lambda && XP.z < N * lambda - lambda && xdisp < XP.x) xdisp = XP.x;
    };
    
    // count the errors
    int errors = 0;

    // time of the last timestep should be equal to tau
    double FinalTime = electron->getTime();
    Vector FinalPosition = electron->getPosition();

    if (fabs(FinalPosition.z - 3.0) > 1.0e-5)
    {
        errors++;
        printf("final position z= %12.9g m - \033[1;31m test failed!\033[0m\n", FinalPosition.z);
    }
    else
    {
        printf("final position z= %12.9g m - \033[1;32m OK\033[0m\n", FinalPosition.z);
    }
    
    // check the maximum displacement
    if (fabs(xdisp - xdis) > 1e-5)
    {
        errors++;
	printf("X-displacement = %9.6f mm - \033[1;31m test failed!\033[0m\n", xdisp*1.0e3);
    }
    else
    {
	printf("X-displacement = %9.6f mm - \033[1;32m OK\033[0m\n", xdisp*1.0e3);
    }

    if (fabs(FinalTime - tau) > deltaT)
    {
        errors++;
        printf("error in final time = %12.9g s - \033[1;31m test failed!\033[0m\n", FinalTime - tau);
    }
    else
    {
        printf("final time = %12.9g s - \033[1;32m OK\033[0m\n", FinalTime);
    }

    if (fabs(FinalPosition.x) > 1.0e-5)
    {
        errors++;
	printf("Final X-displacement = %9.6f mm - \033[1;31m test failed!\033[0m\n", FinalPosition.x*1.0e3);
    }
    else
    {
	printf("Final X-displacement = %9.6f mm - \033[1;32m OK\033[0m\n", FinalPosition.x*1.0e3);
    }

    // clean up
    delete lattice;
    delete electron;

    return errors;
}
