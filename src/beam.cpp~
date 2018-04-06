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

#include "beam.h"
#include "bunch.h"
#include "particle.h"
#include "global.h"
#include <iostream>

Beam::Beam()
{
    NOB = 0;
}

Beam::~Beam()
{
    for(int i=0; i<NOB; i++)
	delete B[i];
}

void Beam::Add(Bunch *bunch)
{
    NOB++;
    B.push_back(bunch);
}

int Beam::getNOB()
{
    return NOB;
}

double Beam::getTotalCharge()
{
    double charge = 0.0;
    for(int i=0; i<NOB; i++)
    {
	charge += B[i]->getTotalCharge();
    }
    return charge;
}

void Beam::InitVay(double tstep,
		GeneralField *field)
{
    dt = tstep;
    for(int i=0; i<NOB; i++)
    {
	B[i]->InitVay(tstep, field);
    }
}

void Beam::StepVay(GeneralField *field)
{
    for(int i=0; i<NOB; i++)
    {
	B[i]->StepVay(field);
    }
}

//! @todo not yet coded - how to store several bunches hierarchically ?
int Beam::WriteWatchPointHDF5(const char *filename)
{
    return 0;
}

void Beam::integrateFieldTrace(
    Vector ObservationPoint,
    double t0,
    double dt,
    int nots,
    std::vector<ElMagField> *ObservationField)
{
    // just do the summation over all te bunches
    for(int i=0; i<NOB; i++)
    {
	B[i]->integrateFieldTrace(ObservationPoint,t0,dt,nots,ObservationField);
    }
}
