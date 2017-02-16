/*=========================================================================
 * 
 *  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers
 * 
 *  Copyright (c) 2017 U. Lehnert
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * =========================================================================*/

#ifndef HOMOGENEOUSMAGNET_H
#define HOMOGENEOUSMAGNET_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "externalfield.h"
#include "vector.h"

using namespace std;

//a class for a homogenous magnetic field 
//the field would extend everywhere

class HomogeneousMagnet : public ExternalField
{

  public:

    // constructor
    HomogeneousMagnet(Vector B);        // field vector B[T]

    Vector getB0();			// report the field vector
    
  private:

    Vector ElementLocalEField(double t, Vector X);

    Vector ElementLocalBField(double t, Vector X);

  private:
      
    Vector B0;				// the constant field vector

};

#endif
