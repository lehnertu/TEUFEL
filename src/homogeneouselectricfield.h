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

#ifndef HOMOGENEOUSELECTRICFIELD_H
#define HOMOGENEOUSELECTRICFIELD_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "externalfield.h"
#include "vector.h"

using namespace std;

//a class for a homogenous electric field 
//the field would extend everywhere

class HomogeneousElectricField : public ExternalField
{

  public:

    // constructor
    HomogeneousElectricField(Vector E);        // field vector E[V/m]

    Vector getE0();			// report the field vector
    
  private:

    Vector ElementLocalEField(double t, Vector X);

    Vector ElementLocalBField(double t, Vector X);

  private:
      
    Vector E0;				// the constant field vector

};

#endif
