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

#include "csr.h"
#include "global.h"

#include <iostream>
#include <math.h>
#include "hdf5.h"

CSR::CSR()
{
    is_initialized = false;
    N_slices = 0;
    createOutput = false;
}

CSR::CSR(
    const pugi::xml_node node,
    InputParser *parser )
{
    // no source definet yet
    is_initialized = false;
    pugi::xml_attribute att = node.attribute("N_slices");
    if (!att)
        throw(IOexception("InputParser::CSR_2D - attribute N_slices not found."));
    N_slices = parser->parseInt(att);
    // define file output if requested
    pugi::xml_node lognode = node.child("log");
    if (lognode)
    {
        pugi::xml_attribute fn = lognode.attribute("file");
        if (!fn) throw(IOexception("InputParser::CSR - <log> filename for log not found."));
        FileName = fn.as_string();
        createOutput = true;
    } else {
        createOutput = false;
    }
    if (teufel::rank==0)
    {
        std::cout << "CSR interaction" << std::endl;
        std::cout << "  " << N_slices << " slices" << std::endl;
    }
}

void CSR::init()
{
    if (is_initialized)
        throw(IOexception("error - CSR::init() called twice."));
    is_initialized = true;
}

CSR::~CSR()
{
    // free field map memory
    // if (is_initialized)
    //    delete[] interaction_field;
    // free output storage memory
}

void CSR::update(Beam *beam, double tracking_time)
{
    double update_time = tracking_time;
    if (teufel::rank==0)
    {
        std::cout << "CSR::update() at tracking time " << update_time << " s" << std::endl;
    }
    // get the particle coordinates
    int NOP = beam->getNOP();
    double *buffer = new double[NOP*6];
    beam->bufferCoordinates(buffer, NOP);
    Vector *position = new Vector[NOP];
    Vector *momentum = new Vector[NOP];
    double *bp = buffer;
    for(int i=0; i<NOP; i++)
    {
        position[i].x = *bp++;
        position[i].y = *bp++;
        position[i].z = *bp++;
        momentum[i].x = *bp++;
        momentum[i].y = *bp++;
        momentum[i].z = *bp++;
    };
    
    delete[] buffer;
    delete[] position;
    delete[] momentum;
}

ElMagField CSR::Field(double t, Vector X)
{
    // TODO: returning zero fields for now
    return ElMagField();
}

void CSR::write_output()
{
    if (createOutput)
    {
        cout << "CSR : writing slice data to " << FileName << endl;

        // TODO: write actual data

        // no errors have occured if we made it 'til here
        cout << "writing HDF5 done." << endl;
    }
}

