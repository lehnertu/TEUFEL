/*=========================================================================

  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers

  Module:    Radiation test case

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
    \brief Radiation at the centre of a loop test case

    @author Ulf Lehnert
    @author Vipul Joshi
    @date 14.2.2017
    @file teufel.loop.cpp
    
    This test case tracks a single electron in a magnetic field and computes the radiation field generated at the centre of the loop of radius a. 
    The electrons moves in the x-z plane and the radiation is observed at its centre. The angular frequency is given by \f$\omega = q B_0 / (\gamma m_e)\f$
    It is finally compared to the following analytical solution for observing point at r and retarded time \f$t_r\f$; 
    \f[
	E(r,t_r) = \frac{q}{4 \pi \epsilon_0 a^2 c^2} \;[\;((a^2 w^2-c^2)cos(\omega t_r) + a c \omega \sin(\omega t_r))\,i\;+\;((a^2 \omega^2 -c^2) \sin(\omega t_r) - a c \omega \cos(\omega t_r))\;k\;]
    \f]
    \f[
	B(r,t_r) = \frac{q \omega}{4 \pi \epsilon_0 a^2 c^2}\; k
    \f]
	
    The program computes the magnetic field at the centre of the loop by assuming the motion of electron to be current flowing in the circle. 
    The solution for this case must reduce to the field generated by a current carrying loop, whose solution is readily given as 
    \f[
	B = \mu_0\;\frac{I}{2 a}
    \f]
    \f[
	I = \frac{\omega q}{2\pi}
    \f]
    
    The following is thus checked for:
    @li the magnetic field is correctly computed.
    @li the kinetic energy has not changed
    @li solution can be reduced to that of current carrying loop
    
    Using the Vay algorithm, the particle will be tracked for 10ns and changes in its path will be noticed.
    
    @return The number of errors encountered in the above list of checks is reported.
    
    @todo A rather large number of time steps is required in order to obtain a
    reasonable accuracy of the compued field. Most likely, the interpolation of the trajectory
    data needs some improvement.
    
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "particle.h"
#include "fields.h"
#include "homogeneousmagnet.h"

int NOTS = 5000;                // number of time steps

int main ()
{
    printf("\nTEUFEL - Radiation at the centre of a circular loop\n");
    Vector B=Vector(0.0, 0.5, 0.0);
    HomogeneousMagnet *mag = new HomogeneousMagnet(B);
    printf("B0 =  %10.3g T\n",B.norm());
    double gamma = 10.0;
    double betagamma= sqrt(gamma*gamma-1.0);
    printf("gamma =  %12.5g\n",gamma);
    printf("beta gamma =  %12.5g\n",betagamma);
    printf("c*p =  %12.5g MeV\n",1e-6*mecsquared*betagamma);
    // double AngFreq = ElementaryCharge*B.norm()/(gamma*9.1e-31);
    double AngFreq = B.norm()/(gamma*mecsquared)*SpeedOfLight*SpeedOfLight;
    printf("Angular Frequency =  %10.3g 1/s\n",AngFreq);
    double Radius= betagamma*0.511e6/(B.norm()*SpeedOfLight);
    printf("Radius R =  %10.3g m\n",Radius);
    double tau = 2*Pi/AngFreq;
    printf("round trip time =  %10.6g s\n",tau);
    double charge = 1.0e-9;
    double ne = charge/ElementaryCharge;
    printf("charge =  %10.3g C   -    %10.3g electrons\n",charge,ne);
    printf("current I =  %10.3g A\n",charge/tau);
    //compute the radiation fields:
    double RadBField=MuNull*charge/tau/(2.0*Radius);
    printf("Expected Radiation Magnetic Field: %12.5g T\n",RadBField);

    // a simple lattice with just the E and B fields
    Lattice *lattice = new Lattice;
    lattice->addElement(mag);
    // one single charge of 0.16 nC
    ChargedParticle *electron = new ChargedParticle(ne,ne);

    // initial position at the origin
    Vector X0 = Vector(0.0, 0.0, 0.0);
    // initial momentum of the particle
    Vector P0 = Vector(0.0,0.0,betagamma);
    
    // track for the time it takes to complete one circle 
    double deltaT = tau / NOTS;
    // electron->TrackVay(NOTS, deltaT, X0, P0, lattice);
    electron->InitVay(0.0, X0, P0, deltaT, lattice);
    for (int i=0; i<NOTS; i++)
	electron->StepVay(lattice);
    
    // count the errors
    int errors = 0;
  
    Vector MPos=electron->TrajPoint(NOTS/2);
    printf("half-loop position = (%10.3g, %10.3g %10.3g) m  ", MPos.x, MPos.y, MPos.z);
    if ((MPos-Vector(-2.0*Radius, 0.0, 0.0)).norm()>1.0e-3)
    {
	errors++;
	printf(" - \033[1;31m test failed!\033[0m\n"); 
    } else {
	printf(" - \033[1;32m OK\033[0m\n");
    }
    
    Vector FPos=electron->TrajPoint(NOTS);
    printf("final position = (%10.3g, %10.3g %10.3g) m  ", FPos.x, FPos.y, FPos.z);
    if ((FPos-X0).norm()>1.0e-4)
    {
	errors++;
	printf(" - \033[1;31m test failed!\033[0m\n"); 
    } else {
	printf(" - \033[1;32m OK\033[0m\n");
    }
    
    // look for the momentum changes
    Vector FMom=electron->TrajMomentum(NOTS);
    if (fabs(FMom.norm()-P0.norm()) > 1.0e-3) {
	errors++;
	printf("Final Momentum = %12.5f - \033[1;31m test failed!\033[0m\n", FMom.norm());
    } else {
	printf("Final Momentum = %12.5f - \033[1;32m OK\033[0m\n",FMom.norm());}

    // calculate the radiation fields at the time the electron finishes the loop
    // the radiation field is observed at the centre of the loop at (R,0,0)
    ElMagField Obs = electron->RetardedField(tau,Vector(-Radius,0.0,0.0));
    
    if (fabs(Obs.B().norm()-RadBField)/RadBField > 1.0e-3) {
	errors++;
	printf("Field observed = %12.5g T \033[1;31m test failed!\033[0m\n", Obs.B().norm());
    } else {
	printf("Field observed = %12.5g T- \033[1;32m OK\033[0m\n",Obs.B().norm());
    }

    return errors;
    
}

