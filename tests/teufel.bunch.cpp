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
    
    Bunch *B = new Bunch();
    Distribution *dist = new Distribution(6, NOP);
    
    // clean up
    delete B;
    delete dist;

    return errors;
}
