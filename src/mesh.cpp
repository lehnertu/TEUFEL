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

#include "bunch.h"
#include "beam.h"
#include "global.h"
#include "observer.h"
#include "mesh.h"

#include <iostream>
#include <math.h>
#include "SDDS.h"
#include "hdf5.h"

MeshedScreen::MeshedScreen(std::string filename)
{
    FileName = filename;
}

void MeshedScreen::integrate(Beam *src)
{
}

void MeshedScreen::integrate(Bunch *src)
{
}

void MeshedScreen::integrate(Lattice *src)
{
}

void MeshedScreen::integrate_mp(Beam *src, unsigned int NumCores, unsigned int CoreId)
{
}

void MeshedScreen::integrate_mp(Bunch *src, unsigned int NumCores, unsigned int CoreId)
{
}

void MeshedScreen::integrate_mp(Lattice *src, unsigned int NumCores, unsigned int CoreId)
{
}

unsigned int MeshedScreen::getBufferSize()
{
    return 0;
}

double* MeshedScreen::getBuffer()
{
    double* buffer = new double[getBufferSize()];
    if (buffer==0)
        throw(IOexception("MeshedScreen::getBuffer() - error allocating memory."));
    else
    {
    };
    return buffer;
}

void MeshedScreen::fromBuffer(double *buffer, unsigned int size)
{
}

void MeshedScreen::WriteMeshedField()
{
}

void MeshedScreen::generateOutput()
{
    WriteMeshedField();
}
