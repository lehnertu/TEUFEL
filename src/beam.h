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

#pragma once

#include "fields.h"
#include "bunch.h"
#include "vector.h"
#include <vector>

using namespace std;

/*!
    \class Beam
    \brief Ensemble of bunches
 
    @author Ulf Lehnert
    @date 19.9.2017
    
    This is a container holding a number of bunches.
    The main purpose is to parallelize the tracking of several buches but
    also tracking of several distinct fractions of a bunch is possible.
    The interface almost exactly resembles that of a single bunch.
 */
class Beam
{

public:

private:

    //! Number of bunches in the beam.
    int NOB;

    //! we store references to all bunches
    vector<Bunch*> B;

};
