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
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "particle.h"
#include "fields.h"
#include "wave.h"

int NOTS = 1000;                // number of time steps

int main ()
{

    printf("\nTEUFEL - radiation push testcase\n");

	GaussianWave *wave = new GaussianWave();
	wave->Setup(1.0e-6, complex<double>(3.23e8,0.0), 1.0);
	ElMagField F = wave->Field(0.0, Vector(0.0, 0.0, 0.0));
    Vector ObsE = F.E();
    printf("E =  (%9.6g,%9.6g,%9.6g) V/m\n",ObsE.x,ObsE.y,ObsE.z);
    double I0 = dot(ObsE,ObsE)/2.0*EpsNull*SpeedOfLight;
    printf("I0 =  %9.6g W/cm²\n",1e4*I0);
    double lambda = wave->getWavelength();
    printf("lambda =  %9.6g µm\n",1e6*lambda);
	double a0 = 0.85e-9*1e6*lambda*sqrt(1e4*I0);
    printf("a0=%9.6g\n",a0);
	
    // a simple lattice with just the inpinging wave
    Lattice *lattice = new Lattice;
    lattice->addElement(wave);
    
    // one single electron
    ChargedParticle *electron = new ChargedParticle();
    
    // initial position at the origin
    Vector X0 = Vector(0.0, 0.0, 0.0);
    // initial momentum of the particle
    Vector P0 = Vector(0.0, 0.0, 0.0);
    Vector A0 = Vector(0.0, 0.0, 0.0);
    electron->initTrajectory(0.0, X0, P0, A0);
    
    // Track the particle for 10 cycles of the wave
    double tau = wave->getWavelength() / SpeedOfLight;
    double deltaT = 10*tau / NOTS;
    electron->InitVay(deltaT, lattice);
    // record maximum displacement in x-direction
    double xmin = 0.0;
    double xmax = 0.0;
    double ymin = 0.0;
    double ymax = 0.0;
    double zmin = 0.0;
    double zmax = 0.0;
    for (int i=0; i<NOTS; i++)
    {
    	electron->StepVay(lattice);
    	Vector x = electron->getPosition();
    	if (x.x<xmin) xmin = x.x;
    	if (x.x>xmax) xmax = x.x;
    	if (x.y<ymin) ymin = x.y;
    	if (x.y>ymax) ymax = x.y;
    	if (x.z<zmin) zmin = x.z;
    	if (x.z>zmax) zmax = x.z;
    };
    double t = electron->getTime();
    Vector XP = electron->getPosition();
    printf("x(%9.6g s) =  (%9.6g,%9.6g,%9.6g) m\n",t,XP.x,XP.y,XP.z);
    printf("x:(%9.6g, %9.6g)   y:(%9.6g, %9.6g)   z:(%9.6g, %9.6g)\n",xmin,xmax,ymin,ymax,zmin,zmax);

    // count the errors
    int errors = 0;
    
    // clean up
    delete lattice;
    delete electron;

    return errors;
}
