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
	printf("E(%6.3f)=%6.3f  E(%6.3f)=%6.3f\n",t1,f1.E().x,t2,f2.E().x);
    };
    double PreviousRetardedTime(Vector ObservationPoint) {return time1;};
    double RetardedTime(Vector ObservationPoint) {return time2;};
    ElMagField PreviousRetardedField(Vector ObservationPoint){return f1;};
    ElMagField RetardedField(Vector ObservationPoint){return f2;};
};

int main ()
{

    printf("\nTEUFEL - Bunch::integrateFieldTrace() test case\n");

    int errors = 0;
    
    Bunch *B = new Bunch();
    Special *p = new Special(0.0, 1.0, 1.5, 0.5);
    B->Add(p);
    printf("Bunch NOP=%d\n", B->getNOP());

    int nots=8;
    double t0=0.1;
    double dt=0.2;
    std::vector<ElMagField> F1(nots);
    B->integrateFieldTrace(Vector(0.0,0.0,0.0),t0,dt,nots,&F1);
    for (int i=0; i<nots; i++)
    {
	printf("t=%6.2f : Ex=%9.6f\n",t0+(i+0.5)*dt,F1[i].E().x);
    }
    if (fabs(F1[0].E().x -1.3)>1e-6) errors++;
    if (fabs(F1[1].E().x -1.1)>1e-6) errors++;
    if (fabs(F1[2].E().x -0.9)>1e-6) errors++;
    if (fabs(F1[3].E().x -0.7)>1e-6) errors++;
    if (fabs(F1[4].E().x -0.275)>1e-6) errors++;
    if (fabs(F1[5].E().x -0.0)>1e-6) errors++;
    if (fabs(F1[6].E().x -0.0)>1e-6) errors++;
    if (fabs(F1[7].E().x -0.0)>1e-6) errors++;
    printf("\n");

    t0=-0.1;
    dt=0.2;
    std::vector<ElMagField> F2(nots);
    B->integrateFieldTrace(Vector(0.0,0.0,0.0),t0,dt,nots,&F2);
    for (int i=0; i<nots; i++)
    {
	printf("t=%6.2f : Ex=%9.6f\n",t0+(i+0.5)*dt,F2[i].E().x);
    }
    if (fabs(F2[0].E().x -0.725)>1e-6) errors++;
    if (fabs(F2[1].E().x -1.3)>1e-6) errors++;
    if (fabs(F2[2].E().x -1.1)>1e-6) errors++;
    if (fabs(F2[3].E().x -0.9)>1e-6) errors++;
    if (fabs(F2[4].E().x -0.7)>1e-6) errors++;
    if (fabs(F2[5].E().x -0.275)>1e-6) errors++;
    if (fabs(F2[6].E().x -0.0)>1e-6) errors++;
    if (fabs(F2[7].E().x -0.0)>1e-6) errors++;
    printf("\n");
    
    t0=-0.7;
    dt=0.2;
    std::vector<ElMagField> F3(nots);
    B->integrateFieldTrace(Vector(0.0,0.0,0.0),t0,dt,nots,&F3);
    for (int i=0; i<nots; i++)
    {
	printf("t=%6.2f : Ex=%9.6f\n",t0+(i+0.5)*dt,F3[i].E().x);
    }
    if (fabs(F3[0].E().x -0.0)>1e-6) errors++;
    if (fabs(F3[1].E().x -0.0)>1e-6) errors++;
    if (fabs(F3[2].E().x -0.0)>1e-6) errors++;
    if (fabs(F3[3].E().x -0.725)>1e-6) errors++;
    if (fabs(F3[4].E().x -1.3)>1e-6) errors++;
    if (fabs(F3[5].E().x -1.1)>1e-6) errors++;
    if (fabs(F3[6].E().x -0.9)>1e-6) errors++;
    if (fabs(F3[7].E().x -0.7)>1e-6) errors++;
    printf("\n");
    delete B;

    Bunch *B2 = new Bunch();
    B2->Add(new Special(1.0, 2.0, 1.0, 1.0));
    B2->Add(new Special(2.0, 2.2, 1.0, 1.0));
    B2->Add(new Special(2.2, 2.4, 1.0, 0.0));
    B2->Add(new Special(2.2, 2.4, 0.0, 1.0));
    B2->Add(new Special(2.4, 3.0, 1.0, 1.0));
    
    printf("Bunch NOP=%d\n", B2->getNOP());
    
    nots=8;
    t0=0.0;
    dt=0.5;
    std::vector<ElMagField> F4(nots);
    B2->integrateFieldTrace(Vector(0.0,0.0,0.0),t0,dt,nots,&F4);
    for (int i=0; i<nots; i++)
    {
	printf("t=%6.2f : Ex=%9.6f\n",t0+(i+0.5)*dt,F4[i].E().x);
    }
    if (fabs(F4[0].E().x -0.0)>1e-6) errors++;
    if (fabs(F4[1].E().x -0.0)>1e-6) errors++;
    if (fabs(F4[2].E().x -1.0)>1e-6) errors++;
    if (fabs(F4[3].E().x -1.0)>1e-6) errors++;
    if (fabs(F4[4].E().x -1.0)>1e-6) errors++;
    if (fabs(F4[5].E().x -1.0)>1e-6) errors++;
    if (fabs(F4[6].E().x -0.0)>1e-6) errors++;
    if (fabs(F4[7].E().x -0.0)>1e-6) errors++;
    printf("\n");
    delete B2;
    
    printf("errors %d\n",errors);
    return errors;
}
