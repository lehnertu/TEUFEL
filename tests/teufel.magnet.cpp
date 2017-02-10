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
    @li the particle has arrived back at the origin
    @li the kinetic energy has not changed
    
    @return The number of errors encountered in the above list of checks is reported.
    
    @author Ulf Lehnert
    @date 10.2.2017
    @file teufel.magnet.cpp
 
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "particle.h"
#include "externalfield.h"

int NOTS = 1000;                // number of time steps

class HomogeneousMagnet : public ExternalField
{

  public:

    // constructor
    HomogeneousMagnet(double B)                 // peak field [T]
    { BPeak = B; };

  private:

    Vector ElementLocalEField(double t, Vector X)
    { return Vector(0.0 ,0.0 ,0.0); };

    Vector ElementLocalBField(double t, Vector X)
    { return Vector(0.0 ,BPeak ,0.0); };

  public:

    double  BPeak;                              // peak field [T]
};


int main ()
{

    printf("\nTEUFEL - homogeneous magnet testcase\n");

    double B = 0.1;
    HomogeneousMagnet *mag = new HomogeneousMagnet(B);
    printf("B =  %9.6g T\n",mag->BPeak);
    double gamma = 10.0;
    double beta = sqrt(1.0-1.0/(gamma*gamma));
    double betagamma= sqrt(gamma*gamma-1.0);
    printf("beta =  %12.9g\n",beta);
    printf("gamma =  %12.9g\n",gamma);
    printf("c*p =  %12.9g MeV\n",1e-6*mecsquared*betagamma);
    double Rgyro = betagamma * mecsquared / SpeedOfLight / B;
    printf("R =  %12.9g m\n",Rgyro);
    double Omega = beta*SpeedOfLight/Rgyro;
    printf("Omega =  %12.9g/s\n",Omega);
    double tau = 2*Pi*Rgyro/(beta*SpeedOfLight);
    printf("tau =  %12.9g s\n",tau);
    double deltaT = tau/NOTS;
    printf("\n");

    Lattice *lattice = new Lattice;
    lattice->addElement(mag);

    ChargedParticle *electron = new ChargedParticle();

    Vector X0 = Vector(0.0, 0.0, 0.0);
    // initial momentum of the particle
    Vector P0 = Vector(betagamma, 0.0, 0.0);
    electron->TrackVay(NOTS, deltaT, X0, P0, lattice);

    Vector MiddlePosition = electron->TrajPoint(NOTS/2);
    printf("middle position =  (%12.9f, %12.9f, %12.9f) m\n", MiddlePosition.x, MiddlePosition.y, MiddlePosition.z);
    double FinalTime = electron->TrajTime(NOTS);
    Vector FinalPosition = electron->TrajPoint(NOTS);
    Vector FinalMomentum = electron->TrajMomentum(NOTS);
    printf("final time =  %12.9g s\n",FinalTime);
    printf("final position =  (%12.9f, %12.9f, %12.9f) m\n", FinalPosition.x, FinalPosition.y, FinalPosition.z);
    printf("final momentum =  (%12.9f, %12.9f, %12.9f)\n", FinalMomentum.x, FinalMomentum.y, FinalMomentum.z);

}
