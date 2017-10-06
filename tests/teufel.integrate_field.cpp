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
    \brief Test case for handling of Bunch::integrateFieldTrace()

    @author Ulf Lehnert
    @date 6.10.2017
    @file teufel.integrate_field.cpp
    
    We derive a Special particle class where we can pre-define which
    values are returned for the generated fields. Then we query
    a Bunch containig just one such particle for the generated field trace
    and check the known values.
    
    @return The number of errors encountered.
    
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "particle.h"
#include "bunch.h"

class Special : public ChargedParticle
{
public:
    double time1, time2;
    ElMagField f1, f2;
    Special(double t1, double t2, double Ex1, double Ex2)
    {
	time1=t1; time2=t2;
	f1=ElMagField(Vector(Ex1,0.0,0.0),Vector(0.0,0.0,0.0));
	f2=ElMagField(Vector(Ex2,0.0,0.0),Vector(0.0,0.0,0.0));
    };
    double PreviousRetardedTime(Vector ObservationPoint) {return time1;};
    double RetardedTime(Vector ObservationPoint) {return time2;};
    ElMagField PreviousRetardedField(Vector ObservationPoint) {return f1;};
    ElMagField RetardedField(Vector ObservationPoint) {return f2;};
};

int main ()
{

    printf("\nTEUFEL - Bunch::integrateFieldTrace() test case\n");

    int errors = 0;
    
    Bunch *B = new Bunch();
    Special *p = new Special(0.0, 1.0, 1.5, 0.5);
    B->Add(p);

    int nots=11;
    double t0=0.1;
    double dt=0.1;
    std::vector<ElMagField> ObservationField(nots);
    B->integrateFieldTrace(Vector(0.0,0.0,0.0),t0,dt,nots,&ObservationField);
    for (int i=0; i<nots; i++)
    {
	printf("t=%6.2f : Ex=%9.6f\n",t0+i*dt,ObservationField[i].E().x);
    }

    delete B;

    return errors;
}
