/*=========================================================================

  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers

  Module:    Electrostatic test case

  Copyright (c) 2019 U. Lehnert

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
    \brief Electrostatic test case

    @author Ulf Lehnert
    @date 18.2.2019
    @file teufel.electrostatic.cpp
    
    This test case tracks puts a single charge at rest into observed space.
    The observed field should be the classical electrostatic field of a point charge.

    The the charge is given a uniform momentum \f$\beta\gamma \approx 10\f$.
    Now the fields should be deformed due to the observation in lab space.
    In addition magnetic fields arise. The fields are checked near the
    starting point of the trajectory testing the back-extrapolation algorithm
    and at a later point where the field is emitted during the motion.
    
    The electrostatic field in a distance \f$r\f$ from a point charge \f$Q\f$
    is \f$E=\frac{Q}{4\pi\epsilon_0 r^2}\f$. In 1 mm distance from a 1pC charge
    this yields E = 8987.55 V/m.
    
    @return The number of errors encountered in the above list of checks is reported.
    
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "particle.h"
#include "fields.h"
#include <iostream>
#include <fstream>

int NOTS = 2000;                // number of time steps

int main ()
{
    printf("\nTEUFEL - Electrostatic\n");

    // a simple lattice with zero E and B fields
    Vector E, B;
    HomogeneousField *lattice = new HomogeneousField(E,B);
    // one single charge of -1 pC
    ChargedParticle *particle = new ChargedParticle(-1.0e-12/ElementaryCharge,1.0);

    // initial position at the origin
    Vector X0 = Vector(0.0, 0.0, 0.0);
    Vector P0 = Vector(0.0, 0.0, 0.0);
    Vector A0 = Vector(0.0, 0.0, 0.0);
    particle->initTrajectory(0.0, X0, P0, A0);
    
    // count the errors
    int errors = 0;

    // compute the reference field value
    // 1mm distance from a 1pC charge
    double Eref = 1e-12 / (4.0*Pi*EpsNull*1.0e-3*1.0e-3);
    
    // check the fiels of the particle at rest
    ElMagField EB;
    EB = particle->RetardedField(0.0, Vector(0.001,0.0,0.0));
    Vector ObsE = EB.E();
    printf("E =  (%9.6g,%9.6g,%9.6g) V/m",ObsE.x,ObsE.y,ObsE.z);
    if ( (Vector(-1.0,0.0,0.0)*Eref - ObsE).norm() < 1.0e-3 ) {
	    printf(" - \033[1;32m OK\033[0m\n");
    } else {
    	errors++;
    	printf(" - \033[1;31m test failed!\033[0m\n");
    }
    EB = particle->RetardedField(0.0, Vector(-0.001,0.0,0.0));
    ObsE = EB.E();
    printf("E =  (%9.6g,%9.6g,%9.6g) V/m",ObsE.x,ObsE.y,ObsE.z);
    if ( (Vector(1.0,0.0,0.0)*Eref - ObsE).norm() < 1.0e-3 ) {
	    printf(" - \033[1;32m OK\033[0m\n");
    } else {
    	errors++;
    	printf(" - \033[1;31m test failed!\033[0m\n");
    }
    EB = particle->RetardedField(0.0, Vector(0.0,0.0,0.001));
    ObsE = EB.E();
    printf("E =  (%9.6g,%9.6g,%9.6g) V/m",ObsE.x,ObsE.y,ObsE.z);
    if ( (Vector(0.0,0.0,-1.0)*Eref - ObsE).norm()/Eref < 1.0e-6 ) {
	    printf(" - \033[1;32m OK\033[0m\n");
    } else {
    	errors++;
    	printf(" - \033[1;31m test failed!\033[0m\n");
    }
    printf("\n");
    
    // create the trajectory of the particle in motion
    double gamma = 10.0;
    double beta = sqrt(1.0-1.0/(gamma*gamma));
    double betagamma= sqrt(gamma*gamma-1.0);
    // create an arbitrary unit vector and generate particle's momentum
    Vector Dir= Vector(2.0,3.0,5.0);
    // make it easy for debugging only
    // Vector Dir= Vector(0.0,0.0,1.0);
    Dir.normalize();
    Vector Perpendicular = cross(Dir,Vector(1.0,0.0,0.0));
    Perpendicular.normalize();
    printf("beta =  %12.9g\n",beta);
    printf("gamma =  %12.9g\n",gamma);
    printf("c*p =  %12.9g MeV\n",1e-6*mecsquared*betagamma);
    printf("beta*gamma= %12.9g\n",betagamma);
    P0 = Dir*betagamma;
    particle->initTrajectory(0.0, X0, P0, A0);
    // track the particle for a little more than 1m
    double tau=4.0e-9;
    double deltaT = tau/NOTS;
    particle->InitVay(deltaT, lattice);
    for (int i=0; i<NOTS; i++)
	    particle->StepVay(lattice);
    printf("\n");

    Vector ErefV = Perpendicular*Eref*-1.0*gamma;
    printf("Eref =  %9.6g V/m (%9.6g,%9.6g,%9.6g) V/m\n",ErefV.norm(),ErefV.x,ErefV.y,ErefV.z);
    Vector BrefV = cross(Dir*beta*SpeedOfLight,Perpendicular)*-1.0e-12*MuNull/(4*Pi)/(0.001*0.001)*gamma;
    printf("Bref =  %9.6g T (%9.6g,%9.6g,%9.6g) T\n",BrefV.norm(),BrefV.x,BrefV.y,BrefV.z);
    printf("\n");

    // check the field near the origin
    Vector ObsPos = Perpendicular*0.001;
    EB = particle->RetardedField(0.0, ObsPos);
    ObsE = EB.E();
    Vector ObsB = EB.B();
    printf("E =  %9.6g V/m (%9.6g,%9.6g,%9.6g) V/m",ObsE.norm(),ObsE.x,ObsE.y,ObsE.z);
    if ( (ErefV-ObsE).norm()/ErefV.norm() < 1.0e-5 ) {
	    printf(" - \033[1;32m OK\033[0m\n");
    } else {
    	errors++;
    	printf(" - \033[1;31m test failed!\033[0m\n");
    }
    printf("B =  %9.6g T (%9.6g,%9.6g,%9.6g) T",ObsB.norm(),ObsB.x,ObsB.y,ObsB.z);
    if ( (BrefV-ObsB).norm()/BrefV.norm() < 1.0e-5 ) {
	    printf(" - \033[1;32m OK\033[0m\n");
    } else {
    	errors++;
    	printf(" - \033[1;31m test failed!\033[0m\n");
    }
    // list the transverse electric field for 0.01mm of trajectory
    printf("trace : ");
    double ETransListNull[11];
    for (int i=0; i<11; i++)
    {
        double l = -0.005+i*0.001;
        EB = particle->RetardedField(0.0, ObsPos+Dir*l);
        ObsE = EB.E();
        ETransListNull[i]=dot(ObsE,Perpendicular);
        printf("%9.6g  ",ETransListNull[i]);
    };
    printf("\n\n");
        
    // check the field near 1m of trajectory
    double ObsTime = 3.0e-9;
    ObsPos = Dir*ObsTime*beta*SpeedOfLight;
    printf("expected at %9.6g s (%9.6g,%9.6g,%9.6g) m\n",ObsTime,ObsPos.x,ObsPos.y,ObsPos.z);
    particle->CoordinatesAtTime(ObsTime,&X0,&P0);
    printf("Trajectory at %9.6g s (%9.6g,%9.6g,%9.6g) m",ObsTime,X0.x,X0.y,X0.z);
    if ( (ObsPos-X0).norm() < 1.0e-6 ) {
	    printf(" - \033[1;32m OK\033[0m\n");
    } else {
    	errors++;
    	printf(" - \033[1;31m test failed!\033[0m\n");
    }
    printf("\n");
    
    ObsPos += Perpendicular*0.001;
    EB = particle->RetardedField(ObsTime, ObsPos);
    ObsE = EB.E();
    ObsB = EB.B();
    printf("E =  %9.6g V/m (%9.6g,%9.6g,%9.6g) V/m",ObsE.norm(),ObsE.x,ObsE.y,ObsE.z);
    if ( (ErefV-ObsE).norm()/ErefV.norm() < 1.0e-3 ) {
	    printf(" - \033[1;32m OK\033[0m\n");
    } else {
    	errors++;
    	printf(" - \033[1;31m test failed!\033[0m\n");
    }
    printf("B =  %9.6g T (%9.6g,%9.6g,%9.6g) T",ObsB.norm(),ObsB.x,ObsB.y,ObsB.z);
    if ( (BrefV-ObsB).norm()/BrefV.norm() < 1.0e-3 ) {
	    printf(" - \033[1;32m OK\033[0m\n");
    } else {
    	errors++;
    	printf(" - \033[1;31m test failed!\033[0m\n");
    }
    // list the transverse electric field for 0.01mm of trajectory
    printf("trace : ");
    double ETransList[11];
    for (int i=0; i<11; i++)
    {
        double l = -0.005+i*0.001;
        EB = particle->RetardedField(ObsTime, ObsPos+Dir*l);
        ObsE = EB.E();
        ETransList[i]=dot(ObsE,Perpendicular);
        printf("%9.6g  ",ETransList[i]);
    };
    printf("\n\n");

    // check equality of both lists
    for (int i=0; i<11; i++)
        if ((ETransList[i]-ETransListNull[i])/ETransListNull[i] > 1e-5)
        {
            printf("Fields differ at trace point No. %d - \033[1;31m test failed!\033[0m\n",i);
          	errors++;
        };

    // clean up
    delete lattice;
    delete particle;

    return errors;
    
}
