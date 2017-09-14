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

#include "bunch.h"

using namespace std;

/*!
 *    \class Distribution
 *    \brief Main class describing/running the whole simulation.
 * 
 *    @author Ulf Lehnert
 *    @date 14.9.2017
 *    
 *    This class provides the main framework for running a simulation.
 *    All input about the configuration is read from an XML input file.
 */
class Simulation
{

public:
    
    /*! The constructor of the simulation.
     */
    Simulation(const char *filename);
    
    //! Destructor only used to free the memory
    ~Simulation();
    
    //! parse the input file and create all described objects
    void parseInput();
    
    //! run the simulation
    void run();
    
    //! again parse the input file and create all requested output files
    void generateOutput();

private:
    
    //! the input file name
    char input[256];

};

