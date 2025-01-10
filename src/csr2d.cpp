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
    pugi::xml_attribute att = node.attribute("N_long");
    if (!att)
        throw(IOexception("InputParser::CSR_2D - attribute N_long not found."));
    N_long = parser->parseInt(att);
    att = node.attribute("N_trans");
    if (!att)
        throw(IOexception("InputParser::CSR_2D - attribute N_trans not found."));
    N_trans = parser->parseInt(att);
    pugi::xml_node vec = node.child("longitudinal");
    if (!vec)
        throw(IOexception("InputParser::CSR_2D - <longitudinal> not found."));
    else
    {
        double x, y, z;
        x = parser->parseDouble(vec.attribute("x"));
        y = parser->parseDouble(vec.attribute("y"));
        z = parser->parseDouble(vec.attribute("z"));
        e_long = Vector(x,y,z);
        e_long.normalize();
    }
    vec = node.child("transversal");
    if (!vec)
        throw(IOexception("InputParser::CSR_2D - <transversal> not found."));
    else
    {
        double x, y, z;
        x = parser->parseDouble(vec.attribute("x"));
        y = parser->parseDouble(vec.attribute("y"));
        z = parser->parseDouble(vec.attribute("z"));
        e_trans = Vector(x,y,z);
        e_trans -= e_long * dot(e_trans,e_long);
        e_trans.normalize();
    }
    // define file output if requested
    pugi::xml_node lognode = node.child("log");
    if (lognode)
    {
        pugi::xml_attribute fn = lognode.attribute("file");
        if (!fn) throw(IOexception("InputParser::CSR_2D - <log> filename for log not found."));
        FileName = fn.as_string();
        step_Output = 1;
        pugi::xml_attribute st = lognode.attribute("step");
        if (st) step_Output = st.as_int();
        createOutput = true;
    } else {
        createOutput = false;
    }
    if (teufel::rank==0)
    {
        std::cout << "CSR-2D interaction" << std::endl;
        std::cout << "  " << N_long << " x " << N_trans << " grid";
        std::cout << " (" << e_long.x << ", " << e_long.y << ", " << e_long.z <<") x";
        std::cout << " (" << e_trans.x << ", " << e_trans.y << ", " << e_trans.z <<")";
        std::cout << std::endl;
    }
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

