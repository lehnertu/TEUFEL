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

    This test case tracks a single electron in a homogeneous magnetic field \f$\vec B(\vec x)=\vec B_0\f$.
    
    This electron moves on an periodic circular trajectory.
    The cyclotron frequency and trajectory radius are compared to known values.
    \f[
	R_{gyro} = \beta \gamma \frac{m_0 c}{B}
    \f]
    \f[
	\Omega_{c} = \frac{e B}{\gamma m_0}
    \f]

    The program computes the trajectory of the electron starting at the coordinate system origin
    with a velocity perpendicular to the field. The magnitude of the velocity is choosen such that
    it is relativistic with \f$\gamma = 10.0\f$. The trajectory radius should be R=0.169613 m in
    a magnetic field with B=0.1 T. The electron ist tracked for an amount of time
    corresponding to one revolution which is . After that it is checked that :
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

int NOTS = 3000;                // number of time steps

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

    printf("\nTEUFEL Version 0.01 U.Lehnert 10/2014\n");
    printf("homogeneous magnet testcase\n\n");

    double B = 0.1;
    HomogeneousMagnet *mag = new HomogeneousMagnet(B);
    printf("B =  %9.6g T\n",mag->BPeak);
    double gamma = 10.0;
    double betagamma= sqrt(gamma*gamma-1.0);
    printf("gamma =  %9.6g\n",sqrt(1.0+betagamma*betagamma));
    printf("c*p =  %9.6g MeV\n",1e-6*mecsquared*betagamma);
    double Rgyro = betagamma * mecsquared / SpeedOfLight / B;
    printf("R =  %9.6g m\n",Rgyro);

    Lattice *lattice = new Lattice;
    lattice->addElement(mag);

    ChargedParticle *electron = new ChargedParticle();

    Vector X0 = Vector(0.0, 0.0, Rgyro);
    // initial momentum of the particle
    Vector P0 = Vector(betagamma, 0.0, 0.0);
    electron->TrackLF(NOTS, 1.0e-12, X0, P0, lattice);

}
