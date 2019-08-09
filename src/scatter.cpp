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

#include "scatter.h"

#include <iostream>
#include <math.h>

Scatterer::Scatterer()
{
}

Scatterer::~Scatterer()
{
}

void Scatterer::integrate(Beam *src)
{
}

void Scatterer::integrate(Bunch *src)
{
}

void Scatterer::integrate(Lattice *src)
{
}

ElMagField Scatterer::Field(double t, Vector X)
{
    return ElMagField(Vector(0.0,0.0,0.0),Vector(0.0,0.0,0.0));
}

