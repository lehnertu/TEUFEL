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
    tracker = TRACKING_NONE;
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

Bunch* Beam::getBunch(int i)
{
    Bunch *b = 0;
    if (i>=0 && i<NOB) b=B[i];
    return b;
}

int Beam::getNOB()
{
    return NOB;
}

int Beam::getNOP()
{
    int nop = 0;
    for(int i=0; i<NOB; i++)
    {
        nop += B[i]->getNOP();
    }
    return nop;
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

void Beam::setupTracking(GeneralField *field)
{
    switch (tracker)
    {
        case TRACKING_NONE:
            throw(IOexception("Beam::setupTracking - no tracking method provided."));
            break;
        case TRACKING_EULER:
            throw(IOexception("Beam::setupTracking - EULER tracking method not yet implemented."));
            break;
        case TRACKING_VAY:
            InitVay(field);
            break;
        default:
            throw(IOexception("Beam::setupTracking - unknown tracking method."));
    }
}

void Beam::doStep(GeneralField *field)
{
    switch (tracker)
    {
        case TRACKING_NONE:
            throw(IOexception("Beam::doStep - no tracking method provided."));
            break;
        case TRACKING_EULER:
            throw(IOexception("Beam::doStep - EULER tracking method not yet implemented."));
            break;
        case TRACKING_VAY:
            StepVay(field);
            break;
        default:
            throw(IOexception("Beam::doStep - unknown tracking method."));
    }
}

void Beam::InitVay(GeneralField *field)
{
    for(int i=0; i<NOB; i++)
    {
	B[i]->InitVay(dt, field);
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
int Beam::WriteWatchPointHDF5(std::string filename)
{
    int nop = 0;
    return nop;
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
