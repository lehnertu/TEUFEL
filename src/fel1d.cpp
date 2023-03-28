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
    prop = Vector(0.0, 0.0, 1.0);
    e_x = Vector(1.0, 0.0, 0.0);
    e_x = Vector(0.0, 1.0, 0.0);
    createOutput = false;
    setup();
}

FEL_1D::FEL_1D(
    double time_step,
    const pugi::xml_node node,
    InputParser *parser )
{
    dt = time_step;
    pugi::xml_attribute att = node.attribute("N");
    if (!att)
        throw(IOexception("InputParser::FEL1D - attribute N not found."));
    N_field = parser->parseInt(att);
    pugi::xml_node vec = node.child("position");
    if (!vec)
        throw(IOexception("InputParser::FEL1D - <position> not found."));
    else
    {
        double x, y, z;
        x = parser->parseDouble(vec.attribute("x"));
        y = parser->parseDouble(vec.attribute("y"));
        z = parser->parseDouble(vec.attribute("z"));
        head = Vector(x,y,z);
    }
    vec = node.child("prop");
    if (!vec)
        throw(IOexception("InputParser::FEL1D - <prop> not found."));
    else
    {
        double x, y, z;
        x = parser->parseDouble(vec.attribute("x"));
        y = parser->parseDouble(vec.attribute("y"));
        z = parser->parseDouble(vec.attribute("z"));
        prop = Vector(x,y,z);
    }
    prop.normalize();
    vec = node.child("pol");
    if (!vec)
        throw(IOexception("InputParser::FEL1D - <pol> not found."));
    else
    {
        double x, y, z;
        x = parser->parseDouble(vec.attribute("x"));
        y = parser->parseDouble(vec.attribute("y"));
        z = parser->parseDouble(vec.attribute("z"));
        e_x = Vector(x,y,z);
    }
    e_x -= prop*dot(e_x,prop);
    e_x.normalize();
    e_y = cross(prop, e_x);
    // define file output if requested
    pugi::xml_node lognode = node.child("log");
    if (lognode)
    {
        pugi::xml_attribute fn = lognode.attribute("file");
        if (!fn) throw(IOexception("InputParser::FEL1D - <log> filename for log not found."));
        fn_Output = fn.as_string();
        step_Output = 1;
        pugi::xml_attribute st = lognode.attribute("step");
        if (st) step_Output = st.as_int();
        createOutput = true;
    } else {
        createOutput = false;
    }
    setup();
}

void FEL_1D::setup()
{
    // proper orthonormalization of (e_x, e_y, prop) is assumed
    dz = prop * (-dt*SpeedOfLight);
    field_E = std::vector<double>(N_field, 0.0);
    // print properties
    if (teufel::rank==0)
    {
        std::cout << "FEL-1D interaction  N = " << N_field << ",   dz = " << dz.norm() << " m" << std::endl;
        std::cout << "  propagation = (" << prop.x << ", " << prop.y << ", " << prop.z << ")" << std::endl;
        std::cout << "  dz = (" << dz.x << ", " << dz.y << ", " << dz.z << ")" << std::endl;
        std::cout << "  e_E = (" << e_x.x << ", " << e_x.y << ", " << e_x.z << ")" << std::endl;
        std::cout << "  e_b = (" << e_y.x << ", " << e_y.y << ", " << e_y.z << ")" << std::endl;
        std::cout << "  logging to " << fn_Output << " every " << step_Output << " steps." << std::endl;
        std::cout << std::endl;
    }
}

FEL_1D::~FEL_1D()
{
}

ElMagField FEL_1D::Field(double t, Vector X)
{
    return ElMagFieldZero;
}
