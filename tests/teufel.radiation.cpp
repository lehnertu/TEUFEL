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

    @author Ulf Lehnert, Vipul Joshi
    @date 14.2.2017
    @file teufel.radiation.cpp
    
    This test case tracks a single electron in a magnetic field and computes the radiation field generated at the centre of the loop of radius a. 
    The electrons moves in the x-z plane and the radiation is observed at its centre. The angular frequency is given by w=qB/gamma*me
    It is finally compared to the following analytical solution for observing point at r and retarded time tr; 
    
	E{r,tr} = (q/(4*pi*epsilon0*a^2*c^2))*[ ((a^2*w^2-c^2)cos(wtr)+acw sin(wtr))i + (a^2w^2 -c^2)*sin(wtr)-acw cos(wtr))k]
    
    
	B{r,tr}=(q*w/(4*pi*epsilon*a*c^2))k
    
    The program computes the magnetic field at the centre of the loop by assuming the motion of electron to be current flowing in the circle. 
    The solution for this case must reduce to the field generated by a current carrying loop, whose solution is readily given as 
    
	B = mu0 I/2*a
    
    The following is thus chaecked for:
    @li the magnetic field is correctly computed.
    @li the kinetic energy has not changed
    @li solution can be reduced to that of current carrying loop by assuming current to be equal to I= w*q/2*Pi.
    Using the Vay algorithm, the particle will be tracked for 10ns and changes in its path will be noticed.
    
    @return The number of errors encountered in the above list of checks is reported.
    
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "particle.h"
#include "externalfield.h"
#include "homogeneousmagnet.h"
#include <iostream>
#include <fstream>
#include <tuple>
int NOTS = 3000;                // number of time steps
using namespace std;
int main ()
{
    printf("\nTEUFEL - Radiation at the centre of a circular loop\n");
    Vector B=Vector(0.0,0.5,0.0);
    HomogeneousMagnet *mag = new HomogeneousMagnet(B);
    printf("B =  %9.6g T\n",B.norm());
    printf("Bx =  %9.6g T\t",B.x);
    printf("By =  %9.6g T\t",B.y);
    printf("Bz =  %9.6g T\n",B.z);
    double gamma = 30.35823388162631;
    double beta = sqrt(1.0-1.0/(gamma*gamma));
    double betagamma= sqrt(gamma*gamma-1.0);
    printf("beta =  %12.9g\n",beta);
    printf("gamma =  %12.9g\n",gamma);
    printf("c*p =  %12.9g MeV\n",1e-6*mecsquared*betagamma);
    printf("beta*gamma= %12.9g\n",beta*gamma);
    double AngFreq = ElementaryCharge*B.norm()/(gamma*9.1e-31);
    printf("Angular Frequency =  %9.6g rad/sec\n",AngFreq);
    double Radius= betagamma*0.511e6/(B.norm()*SpeedOfLight);
    
    //compute the radiation fields:
    double RadBField=ElementaryCharge*AngFreq/(4*Pi*8.85e-12*Radius*SpeedOfLight*SpeedOfLight);
    cout<<"Expected Radiation Magnetic Field: "<<RadBField<<"T"<<endl;


    // a simple lattice with just the E and B fields
    Lattice *lattice = new Lattice;
    lattice->addElement(mag);
    // one single electron
    ChargedParticle *electron = new ChargedParticle();

    // initial position at the origin
    Vector X0 = Vector(0.0, 0.0, 0.0);
    // initial momentum of the particle
    Vector P0 = Vector(0.0,0.0,betagamma);
    
    // approximately the time taken to complete one circle 
    double tau=2.1697404117106876e-09;
    double deltaT = tau/NOTS;
    electron->TrackVay(NOTS, deltaT, X0, P0, lattice);
    // count the errors
    int errors = 0;
  
   //the last trajectory point should coincide with the beginning point X0
    Vector ExpectedPosition=X0;
    printf("Expected Final Position x = %9.6gm; ",ExpectedPosition.x);
    printf("y = %9.6gm; ",ExpectedPosition.y);
    printf("z = %9.6g m\n",ExpectedPosition.z);
    
    // time of the last timestep should be equal to tau
    double FinalTime = electron->TrajTime(NOTS);
    if (fabs(FinalTime-tau) > deltaT) {
	errors++;
	printf("error in final time = %12.9g s - \033[1;31m test failed!\033[0m\n", FinalTime-tau);
    } else {
	printf("final time = %12.9g s - \033[1;32m OK\033[0m\n", FinalTime);
    }

    
    Vector FinalPosition=electron->TrajPoint(NOTS);
    if (fabs(FinalPosition.x-ExpectedPosition.x) >1.0e-6 && fabs(FinalPosition.y-ExpectedPosition.y) >1.0e-6 && fabs(FinalPosition.z-ExpectedPosition.z) >1.0e-6) {
	errors++;
	printf("final x position = %12.9g m - \033[1;31m test failed!\033[0m\n", FinalPosition.x); 
	printf("final y position = %12.9g m - \033[1;31m test failed!\033[0m\n", FinalPosition.y);
	printf("final z position = %12.9g m - \033[1;31m test failed!\033[0m\n", FinalPosition.z);
    } else {
	printf("final x position = %12.9g m - \033[1;32m OK\033[0m\n",FinalPosition.x);
	printf("final y position = %12.9g m - \033[1;32m OK\033[0m\n",FinalPosition.y);
	printf("final z position = %12.9g m - \033[1;32m OK\033[0m\n",FinalPosition.z);
    }
    // look for the momentum changes
    Vector ExpectedMomentum=P0;
    Vector FinalMomentum=electron->TrajMomentum(NOTS);
    if (ExpectedMomentum.norm()-FinalMomentum.norm() > 1.5e-3) {
	errors++;
	printf("Final Momentum = %9.6f - \033[1;31m test failed!\033[0m\n", FinalMomentum.norm());
    } else {
	printf("Final Momentum = %9.6f - \033[1;32m OK\033[0m\n",FinalMomentum.norm());}

    //calculate the radiation fields
    //the radiation fields has to be observed at the centre of the loop. The circle has the centre at (R,0,0)
    tuple<Vector,Vector> M=electron->RetardedEField(2.1697404117106876e-09/4.0,Vector(Radius,0.0,0.0));
    Vector EFIELD=get<0>(M);
    Vector BFIELD=get<1>(M);
    
    if (BFIELD.norm()-RadBField > 1.5e-19) {
	errors++;
	printf("Radiation Emitted = %9.6fe-19 T- \033[1;31m test failed!\033[0m\n", BFIELD.norm()*1e19);
    } else {
	printf("Radiation Emitted = %9.6fe-19 T- \033[1;32m OK\033[0m\n",BFIELD.norm()*1e19);}

    return errors;
    
}
