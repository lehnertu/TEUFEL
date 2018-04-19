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
    \brief Homogeneous magnetic field test case

    @author Ulf Lehnert
    @author Vipul Joshi
    @date 10.2.2017
    @file teufel.magnet.cpp
    
    This test case tracks a single electron in a homogeneous magnetic field.
    The electron moves on an periodic circular trajectory.
    The cyclotron frequency and trajectory radius are compared to known values.
    \f[
	R_{gyro} = \beta \gamma \frac{m_0 c}{e B}
    \f]
    \f[
	\Omega_{c} = \frac{e B}{\gamma m_0}
    \f]

    The program computes the trajectory of the electron starting at the coordinate system origin
    with a velocity perpendicular to the field. The magnitude of the velocity is choosen such that
    it is relativistic with \f$\gamma = 10.0\f$. The trajectory radius should be R=0.169613 m in
    a magnetic field with B=0.1 T. The electron ist tracked for an amount of time
    corresponding to one revolution which is 3.57273 ns. After that it is checked that :
    @li the time is correct
    @li the maximum distance from the origin is twice the trajectory radius
    @li the particle has arrived back at the origin
    @li the kinetic energy has not changed
    
    Using the Vay algorithm, the tracking reaches the required accuracy for 1000 Steps.
    
    @return The number of errors encountered in the above list of checks is reported.
    
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "particle.h"
#include "fields.h"

int NOTS = 500;                // number of time steps

int main ()
{

    printf("\nTEUFEL - homogeneous magnet testcase\n");

    Vector B=Vector(0.033166247903554,0.05,0.08);
    HomogeneousField *mag = new HomogeneousField(Vector(0.0,0.0,0.0),B);
    Vector B0 = mag->Field(0.0,Vector(0.0,0.0,0.0)).B();
    printf("B =  %9.6g T\n",B0.norm());
    printf("Bx = %9.6g T  ",B0.x);
    printf("By = %9.6g T  ",B0.y);
    printf("Bz = %9.6g T\n",B0.z);
    
    double gamma = 10.0;
    // define a second vector which is normal to the field and the initial velocity
    Vector n= Vector(0.08,0.05,0.03);
    //define the velocity vector to be normal to the magnetic field.
    Vector p=cross(n,B);
    double beta = sqrt(1.0-1.0/(gamma*gamma));
    double betagamma= sqrt(gamma*gamma-1.0);
    printf("beta =  %12.9g\n",beta);
    printf("gamma =  %12.9g\n",gamma);
    printf("c*p =  %12.9g MeV\n",1e-6*mecsquared*betagamma);
    double Rgyro = betagamma * mecsquared / SpeedOfLight / B.norm();
    printf("R =  %12.9g m\n",Rgyro);
    double Omega = beta*SpeedOfLight/Rgyro;
    printf("Omega =  %12.9g/s\n",Omega);
    double tau = 2*Pi*Rgyro/(beta*SpeedOfLight);
    printf("tau =  %12.9g s\n",tau);
    printf("\n");

    // a simple lattice with just the homogeneous dipole field
    Lattice *lattice = new Lattice;
    lattice->addElement(mag);
    
    // one single electron
    ChargedParticle *electron = new ChargedParticle();
    
    // initial potion at the origin
    Vector X0 = Vector(0.0, 0.0, 0.0);
    // initial momentum of the particle
    Vector P0 = p*betagamma/p.norm();
    Vector A0 = Vector(0.0, 0.0, 0.0);
    electron->initTrajectory(0.0, X0, P0, A0);
    
    // count the errors
    int errors = 0;
    
    // track the particle for the duration of one half revolution
    double deltaT = 0.5*tau/NOTS;
    electron->InitVay(deltaT, lattice);
    for (int i=0; i<NOTS; i++)
	electron->StepVay(lattice);
    
    // the position after half the time steps is the farthest from origin - twice the radius
    Vector MiddlePosition = electron->getPosition();
    double dist = MiddlePosition.norm();
    if (fabs(0.5*dist-Rgyro) > 1e-5) {
	errors++;
	printf("trajectory radius = %9.6f m - \033[1;31m test failed!\033[0m\n", 0.5*dist);
    } else {
	printf("trajectory radius = %9.6f m - \033[1;32m OK\033[0m\n", 0.5*dist);
    }

    // track the particle for the second half revolution
    for (int i=0; i<NOTS; i++)
	electron->StepVay(lattice);
    
    Vector XF = electron->getPosition();
    printf("final position =  (%9.6f, %9.6f, %9.6f) m\n",XF.x, XF.y, XF.z);
    
    // time of the last timestep should be equal to tau
    double FinalTime = electron->getTime();
    if (1e12*fabs(FinalTime-tau) > 0.01) {
	errors++;
	printf("final time = %12.9g s - \033[1;31m test failed!\033[0m\n", FinalTime);
    } else {
	printf("final time = %12.9g s - \033[1;32m OK\033[0m\n", FinalTime);
    }
    // the position after the last timestep should be at the origin again
    Vector FinalPosition = electron->getPosition();
    if (FinalPosition.norm() > 1e-5) {
	errors++;
	printf("final position error = %6.2f µm - \033[1;31m test failed!\033[0m\n", 1e6*FinalPosition.norm());
    } else {
	printf("final position error = %6.2f µm - \033[1;32m OK\033[0m\n", 1e6*FinalPosition.norm());
    }
    // final momentum should be the same as at start
    Vector FinalMomentum = electron->getMomentum();
    if (fabs(FinalMomentum.norm()-betagamma) > 1e-4) {
	errors++;
	printf("beta gamma = %9.6f - \033[1;31m test failed!\033[0m\n", FinalMomentum.norm());
    } else {
	printf("beta gamma = %9.6f - \033[1;32m OK\033[0m\n", FinalMomentum.norm());
    }	
    printf("\n");

    // clean up
    delete lattice;
    delete electron;

    return errors;
}
