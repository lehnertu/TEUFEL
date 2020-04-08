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
   
    We derive a special particle class from ChrgedParticle. 
    For this class we can custom-define the time and field vectors
    that are returned by RetardedTime() and RetardedField(). Then we query
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
#include "fieldtrace.h"

class Special : public ChargedParticle
{
public:
    Special(double t1, double t2, double Ex1, double Ex2)
    {
        // our virtual trajectory has 2 points
        // this must be set for the loop in integrateFieldTrace() to work correctly
        NP=2;
        // define the two trajectory points (their retarded fields)
	    time1=t1; time2=t2;
	    f1=ElMagField(Vector(Ex1,0.0,0.0),Vector(0.0,0.0,0.0));
	    f2=ElMagField(Vector(Ex2,0.0,0.0),Vector(0.0,0.0,0.0));
	    printf("E(%6.3f)=%6.3f  E(%6.3f)=%6.3f\n",t1,f1.E().x,t2,f2.E().x);
    };
    virtual double RetardedTime(int index, Vector ObservationPoint)
    {
        if (index==0)
            return time1;
        else
            return time2;
    };
    virtual ElMagField RetardedField(int index, Vector ObservationPoint)
    {
        if (index==0)
            return f1;
        else
            return f2;
    };
private:
    double time1, time2;
    ElMagField f1, f2;
};

int main ()
{

    printf("\nTEUFEL - Bunch::integrateFieldTrace() test case\n");

    int errors = 0;
    
    Bunch *B = new Bunch();
    Special *p = new Special(0.0, 1.0, 1.5, 0.5);
    B->Add(p);
    printf("Bunch NOP=%d\n\n", B->getNOP());

    // a FieldTrace object aligned with the limits of the trace
    FieldTrace *ft = new FieldTrace(0.0, 0.2, 8);
    B->integrateFieldTrace(Vector(0.0,0.0,0.0),ft);
    for (std::size_t i=0; i<8; i++)
    {
    	printf("t=%6.2f : Ex=%9.6f\n",ft->get_time(i),ft->get_field(i).E().x);
    }
    if (fabs(ft->get_field((std::size_t)0).E().x - 0.725)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)1).E().x - 1.3)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)2).E().x - 1.1)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)3).E().x - 0.9)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)4).E().x - 0.7)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)5).E().x - 0.275)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)6).E().x - 0.0)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)7).E().x - 0.0)>1e-6) errors++;
    printf("\n");
    delete ft;
    
    // now a FieldTrace object not aligned to the track endpoints
    ft = new FieldTrace(-0.1, 0.2, 8);
    B->integrateFieldTrace(Vector(0.0,0.0,0.0),ft);
    for (std::size_t i=0; i<8; i++)
    {
    	printf("t=%6.2f : Ex=%9.6f\n",ft->get_time(i),ft->get_field(i).E().x);
    }
    if (fabs(ft->get_field((std::size_t)0).E().x - 0.0)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)1).E().x - 1.4)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)2).E().x - 1.2)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)3).E().x - 1.0)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)4).E().x - 0.8)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)5).E().x - 0.6)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)6).E().x - 0.0)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)7).E().x - 0.0)>1e-6) errors++;
    printf("\n");
    delete ft;
    delete B;
    
    // a short field spike
    B = new Bunch();
    p = new Special(0.05, 0.15, 1.5, 0.5);
    B->Add(p);
    printf("Spike NOP=%d\n\n", B->getNOP());

    ft = new FieldTrace(-0.2, 0.2, 4);
    B->integrateFieldTrace(Vector(0.0,0.0,0.0),ft);
    for (std::size_t i=0; i<4; i++)
    {
    	printf("t=%6.2f : Ex=%9.6f\n",ft->get_time(i),ft->get_field(i).E().x);
    }
    if (fabs(ft->get_field((std::size_t)0).E().x - 0.0)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)1).E().x - 0.3125)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)2).E().x - 0.1875)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)3).E().x - 0.0)>1e-6) errors++;
    printf("\n");
    delete ft;

    ft = new FieldTrace(-0.1, 0.2, 4);
    B->integrateFieldTrace(Vector(0.0,0.0,0.0),ft);
    for (std::size_t i=0; i<4; i++)
    {
    	printf("t=%6.2f : Ex=%9.6f\n",ft->get_time(i),ft->get_field(i).E().x);
    }
    if (fabs(ft->get_field((std::size_t)0).E().x - 0.0)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)1).E().x - 0.5)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)2).E().x - 0.0)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)3).E().x - 0.0)>1e-6) errors++;
    printf("\n");
    delete ft;
    delete B;
    
    // an extended field trace simulated with 5 particles partially overlapping
    printf("\n");
    B = new Bunch();
    B->Add(new Special(1.0, 2.0, 1.0, 0.0));
    B->Add(new Special(1.0, 2.0, 0.0, 1.0));
    B->Add(new Special(2.0, 2.2, 1.0, 1.0));
    B->Add(new Special(2.2, 2.4, 1.0, 0.0));
    B->Add(new Special(2.4, 3.0, 0.0, 1.0));
    printf("Bunch NOP=%d\n\n", B->getNOP());
    
    ft = new FieldTrace(0.0, 0.5, 8);
    B->integrateFieldTrace(Vector(0.0,0.0,0.0),ft);
    for (std::size_t i=0; i<8; i++)
    {
    	printf("t=%6.2f : Ex=%9.6f\n",ft->get_time(i),ft->get_field(i).E().x);
    }
    if (fabs(ft->get_field((std::size_t)0).E().x - 0.0)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)1).E().x - 0.0)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)2).E().x - 0.5)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)3).E().x - 1.0)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)4).E().x - 0.9875)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)5).E().x - 0.316667)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)6).E().x - 0.395833)>1e-6) errors++;
    if (fabs(ft->get_field((std::size_t)7).E().x - 0.0)>1e-6) errors++;
    delete ft;
    delete B;

    if (errors>0)
        printf("\033[1;31m errors %d - test failed!\033[0m\n",errors);
    else
    	printf("\033[1;32m errors %d - OK\033[0m\n",errors);
    
    return errors;
}
