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
    \brief Test case for handling of particle bunches

    @author Ulf Lehnert
    @date 17.5.2017
    @file teufel.bunch.cpp
    
    Create a particle bunch, copy it and delete the original bunch.
    Modify the copy and output it to file.
    Checks for memory leaks can be performed using valgrind.
    
    @return The number of errors encountered.
    
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "particle.h"
#include "bunch.h"

int main ()
{

    printf("\nTEUFEL - bunch test case\n");

    // number of particles
    int NOP = 100000;

    // count the errors
    int errors = 0;

    double gamma = 50.0;
    double betagamma = sqrt(gamma * gamma - 1.0);
    Distribution *dist = new Distribution(6, NOP);
    dist->generateGaussian(0, 0.002);	// x gaussian displaced 1mm with sigma=2mm
    dist->generateGaussian(1, 0.001);	// y gaussian with sigma=1mm
    dist->generateGaussian(2, 0.0003);	// z gaussian with sigma=0.3mm
    dist->generateGaussian(3, 0.001);	// px gaussian
    dist->generateGaussian(4, 0.001);	// py gaussian
    dist->generateGaussian(5, 0.001);	// pz gaussian
    Bunch *B = new Bunch(dist, 0.0, Vector(0.001,0.0,0.0), Vector(0.0,0.0,betagamma), -1.0, 1.0);

    for (int i=0; i<6; i++)
    {
	ChargedParticle *P = B->getParticle(i);
	Vector X = P->getPosition();
	printf("(%d): x=%g",i,X.x);
	printf(" z=%g\n",X.z);
	Vector BG = P->getMomentum();
	printf("(%d): bg.x=%g",i,BG.x);
	printf(" bg.z=%g\n",BG.z);
    }
    
    // create a particle dump
    if (0 != B->WriteWatchPointSDDS("teufel_bunch_start.sdds"))
    {
	errors++;
	printf("SDDS write \033[1;31m failed!\033[0m\n");
    }
    else
    {
	printf("SDDS file written - \033[1;32m OK\033[0m\n");
    }

    // clean up
    delete B;
    delete dist;

    return errors;
}
