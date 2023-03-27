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

/*!
 * \brief defines the natural constants
 * 
 * @author Ulf Lehnert
 * @date 1.1.2012
 * @file global.h
 */

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <exception>

#define ElementaryCharge 1.6021766208e-19       // [As]
#define mecsquared 0.5109989461e6               // electron mass [eV]
#define SpeedOfLight 2.99792458e8               // [m/s]
#define PlanckH 6.62607015e-34                  // h [J s]
#define Pi 3.1415926535897932384626433832795029
#define EpsNull 8.854187817e-12                 // vacuum permittivity [As/Vm] 1/c² = e0*µ0
#define MuNull (4*Pi*1e-7)                      // magnetic field constant [Vs/Am]
#define InvRestMass 1.758820025e11              // 1 / m = c² / mc² [m²/s²/eV]

// the electron rest mass in kg can be obtained as
// mecsquared*ElementaryCharge/SpeedOfLight^2

/*! define debug levels:
 *  0: suppress all printouts except for those leading to an immediate abort
 *  1: show errors and summary reports
 *  2: show all warnings
 *  3: print additional debugging information
 */
#define DEBUGLEVEL 1

namespace teufel
{
    /*! The rank denotes which of the many nodes of an MPI cluster the code is running on.
     *  It is necessary to know the tasks this particular node has to accomplish
     *  The MPI rank is kept as a global variable accessible from all modules.
     *  It is set to values different from 0 only by MPI executables.
     *
     *  Must be set to zero for all non-MPI executables.
     *  It is used to suppress unnecessary output from code
     *  executed by many nodes in parallel.
     */
    extern int rank;
}

//! time in seconds since the Epoch
double current_time();

/*! 
 * \class IOexception
 * \brief Class for exceptions returning a user-defined message.
 */
class IOexception: public std::exception
{
public:
    IOexception(const char* message) {m=message;}
    virtual const char* what() const throw() {return m;};
    const char* m;
};

/*! 
 * \class RANGEexception
 * \brief Class for exceptions returning a user-defined message.
 */
class RANGEexception: public std::exception
{
public:
    RANGEexception(const char* message) {m=message;}
    virtual const char* what() const throw() {return m;};
    const char* m;
};

/*! a data struct to hold the definition of a watch point
 *  where particle coordinates will be written during tracking
 */
enum WatchType { hdf5, sdds };
typedef struct {
    int step;
    std::string filename;
    WatchType type;
    } watch_t;

