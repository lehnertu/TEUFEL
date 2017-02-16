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
	x_max{undu} = e*B/(gamma*me*beta*c*ku^2)
    \f]
    \f[
	slippage ~ radiation wavelength
    \f]

    The program computes the trajectory of the electron starting at the coordinate system origin
    with a velocity in z direction. The magnitude of the velocity is choosen such that
    it is relativistic with \f$\gamma = 30.35823\f$. The trajectory's maximum displacment in x should be R=0.35686mm in
    a magnetic field with B=0.532 T.  After tracking the particle, it is checked that :
    @li the maximum displacement is correct
    @li the slippage is very close to the radiation wavelength
    @li the particle has returned to x=0,y=0 and moving along the z-axis
    @li the kinetic energy has not changed
    @li the electron is out of the undulator
    
    Using the Vay algorithm, the tracking reaches the required accuracy for 1000 Steps.
    
    @return The number of errors encountered in the above list of checks is reported.
    
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "particle.h"
#include "externalfield.h"
#include "undulator.h"
#include <iostream>
#include <fstream>

int NOTS = 3000;                // number of time steps

int main ()
{
    printf("\nTEUFEL - Undulator testcase\n");
    double B=0.532;
    double lambda=0.037;
    double N=54;
    double ku=2*Pi/lambda;
    Undulator *Undu = new Undulator(B,lambda,N);
    printf("B =  %9.6g T\n",B);
    printf("Undulator Period = %9.6g T\n ",lambda);
    printf("N = %9.6g T\n ",(double)N);   
    double gamma = 30.35823388162631;
    double beta = sqrt(1.0-1.0/(gamma*gamma));
    double betagamma= sqrt(gamma*gamma-1.0);
    double K=0.934*100*lambda*B;
    double lambdar=(lambda/(2*gamma*gamma))*(1+K*K/2);
    printf("beta =  %12.9g\n",beta);
    printf("gamma =  %12.9g\n",gamma);
    printf("c*p =  %12.9g MeV\n",1e-6*mecsquared*betagamma);
    double xdis = ElementaryCharge*B*SpeedOfLight/(gamma*beta*9.1e-31*SpeedOfLight*SpeedOfLight*ku*ku);
    printf("X-Displacement =  %12.9g m\n",xdis);
    printf("Radiation Wavelength =  %12.9g/s\n",lambdar);
    double avgvz=(1-(1/(2*gamma*gamma))*(1+K*K/2))*SpeedOfLight;


    // a simple lattice with just the Undulator Field
    Lattice *lattice = new Lattice;
    lattice->addElement(Undu);

    // one single electron
    ChargedParticle *electron = new ChargedParticle();

    // initial position at the origin
    Vector X0 = Vector(0.0, 0.0, -lambda/2);
    // initial momentum of the particle
    Vector P0 = Vector(0.0,0.0,betagamma);
    
    // track the particle for the duration of single pass
    double tau=06.8e-9;
    double deltaT = tau/NOTS;
    electron->TrackVay(NOTS, deltaT, X0, P0, lattice);

    // count the errors
    int errors = 0;
    
    // time of the last timestep should be equal to tau
    double FinalTime = electron->TrajTime(NOTS);
    Vector FinalPosition=electron->TrajPoint(NOTS);
    double xdisp=0; //maximum displacement in + direction
    if (FinalPosition.z < N*lambda) {
	errors++;
	printf("final position = %12.9g s - \033[1;31m test failed!\033[0m\n", FinalPosition.z);
    } else {
	printf("final position = %12.9g s - \033[1;32m OK\033[0m\n", FinalPosition.z);
    }
    // look for the maximum displacement
    for (int i=0;i<NOTS;i++)
	{
		Vector XP=electron->TrajPoint(i);
		if(XP.z>lambda && XP.z<N*lambda-lambda && xdisp<XP.x) //check for the maximum displacement after the endpieces
		{
			xdisp=XP.x;
		}
	}
    if (fabs(xdisp-xdis) > 1e-5) {
	errors++;
	printf("X-displacement = %9.6f m - \033[1;31m test failed!\033[0m\n", xdisp);
    } else {
	printf("X-displacement = %9.6f m - \033[1;32m OK\033[0m\n", xdisp);
    }
    if (fabs(FinalTime-tau) > deltaT) {
	errors++;
	printf("error in final time = %12.9g s - \033[1;31m test failed!\033[0m\n", FinalTime-tau);
    } else {
	printf("final time = %12.9g s - \033[1;32m OK\033[0m\n", FinalTime);
    }

    if (fabs(electron->TrajPoint(NOTS).x) >1.5e-6 ) {
	errors++;
	printf("Final X-displacement = %9.6f m - \033[1;31m test failed!\033[0m\n",electron->TrajPoint(NOTS).x);
    } else {
	printf("Final X-displacement = %9.6f m - \033[1;32m OK\033[0m\n", electron->TrajPoint(NOTS).x);
    }

    //compute the total slippage between the photon and the electron inside the undulator
    double ds=0; //path length of the electron
    //double ds1=0; //path length of the photon
    //look for the zero crossover after the first endpiece and before the first endpiece.
    //compute average velocity of the electron inside the undulator
    //betaz actually represents the velocity of the electron in the z direction.do not confuse with the beta.
    int j=0;
    int k=0;
    for (int i=0;i<NOTS;i++)
	{
	   Vector XP=electron->TrajPoint(i);
	   Vector XP1=electron->TrajPoint(i+1);
	   if (XP.z<lambda && XP.x*XP1.x <0 )
		{
		   j= i+1;
		}
	   if (XP.z>lambda && XP.z<2*lambda && XP.x*XP1.x<0  )
		{
		  k=i;
		}
	}
double dz=0;  
    for (int i=j;i<k;i++)
    {
        double x1=electron->TrajPoint(i+1).x;
	double x2=electron->TrajPoint(i).x;
	double z1=electron->TrajPoint(i+1).z;
	double z2=electron->TrajPoint(i).z;
	double Vz=(electron->TrajMomentum(i).z);
	double Vz1=(electron->TrajMomentum(i+1).z);
	double betaz=(sqrt(Vz*Vz/(1+Vz*Vz))+sqrt(Vz1*Vz1/(1+Vz1*Vz1)))*SpeedOfLight/2.0;
	
	ds=ds+sqrt(1+((x2-x1)/(z2-z1))*((x2-x1)/(z2-z1)))*(betaz*deltaT);
	dz=dz+(z2-z1);
    }
    
    double slippage=SpeedOfLight*ds/avgvz-ds;
    if (fabs(slippage-lambdar) >1.5e-6 ) {
	errors++;
	printf("Slippage after one period= %9.6f m - \033[1;31m test failed!\033[0m\n",slippage);
    } else {
	printf("Slippage after one period= %9.6f m - \033[1;32m OK\033[0m\n", slippage);
    }
    
    
    return errors;
}
