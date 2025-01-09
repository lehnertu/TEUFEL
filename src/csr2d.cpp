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

#include "csr2d.h"
#include "global.h"

#include <iostream>
#include <math.h>
#include "hdf5.h"

CSR_2D::CSR_2D()
{
    N_long = 0;
    N_trans = 0;
    e_long = Vector(0.0, 0.0, 1.0);
    e_trans = Vector(1.0, 0.0, 0.0);
    N_stored = 0;
    createOutput = false;
    step_Output = 1;
}

CSR_2D::CSR_2D(
    const pugi::xml_node node,
    InputParser *parser )
{
}

void CSR_2D::init(Beam *beam)
{
    source = beam;
    // TODO: sanity checks, timing setup, ...
}

void CSR_2D::update(double tracking_time, double tracking_time_step)
{
    // TODO
}

ElMagField CSR_2D::Field(double t, Vector X)
{
    // TODO: return zero fields for now
    return ElMagField();
}

void CSR_2D::write_output()
{
    if (createOutput)
    {
    }
}
    
CSR_2D::~CSR_2D()
{
    // TODO: free storage
}

