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

#include "fel1d.h"
#include "global.h"

#include <math.h>

FEL_1D::FEL_1D(
    double time_step,
    int number_steps )
{
    N_field = number_steps;
    dt = time_step;
}

FEL_1D::FEL_1D(
    double time_step,
    const pugi::xml_node node,
    InputParser *parser )
{
    dt = time_step;
}

FEL_1D::~FEL_1D()
{
}

ElMagField FEL_1D::Field(double t, Vector X)
{
    return ElMagFieldZero;
}
