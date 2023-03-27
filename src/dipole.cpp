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

#include "dipole.h"
#include "global.h"

#include <math.h>

HardEdgeDipole::HardEdgeDipole() :
    LocalizedField()
{
    B_val = Vector(0.0, 0.0, 0.0);
    p_1 = Vector(0.0, 0.0, 0.0);
    n_1 = Vector(1.0, 0.0, 0.0);
    p_2 = Vector(0.0, 0.0, 0.0);
    n_2 = Vector(-1.0, 0.0, 0.0);
}

HardEdgeDipole::HardEdgeDipole(Vector B, Vector p1, Vector n1, Vector p2, Vector n2) :
    LocalizedField()
{
    B_val = B;
    p_1 = p1;
    n_1 = n1;
    p_2 = p2;
    n_2 = n2;
}

HardEdgeDipole::HardEdgeDipole(const pugi::xml_node node, InputParser *parser) :
    LocalizedField()
{
    parser->parseCalcChildren(node);
    pugi::xml_node field = node.child("B");
    if (!field)
        throw(IOexception("InputParser::HardEdgeDipole - <B> not found."));
    else
    {
        double x, y, z;
        x = parser->parseDouble(field.attribute("x"));
        y = parser->parseDouble(field.attribute("y"));
        z = parser->parseDouble(field.attribute("z"));
        B_val = Vector(x,y,z);
    }
    pugi::xml_node pos1 = node.child("p1");
    if (!pos1)
        throw(IOexception("InputParser::HardEdgeDipole - <p1> not found."));
    else
    {
        double x, y, z;
        x = parser->parseDouble(pos1.attribute("x"));
        y = parser->parseDouble(pos1.attribute("y"));
        z = parser->parseDouble(pos1.attribute("z"));
        p_1 = Vector(x,y,z);
    }
    pugi::xml_node norm1 = node.child("n1");
    if (!norm1)
        throw(IOexception("InputParser::HardEdgeDipole - <n1> not found."));
    else
    {
        double x, y, z;
        x = parser->parseDouble(norm1.attribute("x"));
        y = parser->parseDouble(norm1.attribute("y"));
        z = parser->parseDouble(norm1.attribute("z"));
        n_1 = Vector(x,y,z);
    }
    pugi::xml_node pos2 = node.child("p2");
    if (!pos2)
        throw(IOexception("InputParser::HardEdgeDipole - <p2> not found."));
    else
    {
        double x, y, z;
        x = parser->parseDouble(pos2.attribute("x"));
        y = parser->parseDouble(pos2.attribute("y"));
        z = parser->parseDouble(pos2.attribute("z"));
        p_2 = Vector(x,y,z);
    }
    pugi::xml_node norm2 = node.child("n2");
    if (!norm2)
        throw(IOexception("InputParser::HardEdgeDipole - <n2> not found."));
    else
    {
        double x, y, z;
        x = parser->parseDouble(norm2.attribute("x"));
        y = parser->parseDouble(norm2.attribute("y"));
        z = parser->parseDouble(norm2.attribute("z"));
        n_2 = Vector(x,y,z);
    }
}

ElMagField HardEdgeDipole::LocalField(double t, Vector X)
{
    Vector E = Vector(0.0, 0.0, 0.0);
    Vector B = Vector(0.0, 0.0, 0.0);
    // check distances from the limiting planes
    double d1 = dot(X-p_1,n_1);
    double d2 = dot(X-p_2,n_2);
    // non-zero part of the field
    if ( (d1<0.0) && (d2<0.0) ) B = B_val;
    return ElMagField(E, B);
}




SoftEdgeDipole::SoftEdgeDipole() :
    LocalizedField()
{
    B_val = Vector(0.0, 0.0, 0.0);
    p_1 = Vector(0.0, 0.0, 0.0);
    n_1 = Vector(1.0, 0.0, 0.0);
    t_1 = 0.01;
    p_2 = Vector(0.0, 0.0, 0.0);
    n_2 = Vector(-1.0, 0.0, 0.0);
    t_2 = 0.01;
}

SoftEdgeDipole::SoftEdgeDipole(Vector B, Vector p1, Vector n1, double t1, Vector p2, Vector n2, double t2) :
    LocalizedField()
{
    B_val = B;
    p_1 = p1;
    n_1 = n1;
    t_1 = t1;
    p_2 = p2;
    n_2 = n2;
    t_2 = t2;
}

SoftEdgeDipole::SoftEdgeDipole(const pugi::xml_node node, InputParser *parser) :
    LocalizedField()
{
    parser->parseCalcChildren(node);
    pugi::xml_node field = node.child("B");
    if (!field)
        throw(IOexception("InputParser::SoftEdgeDipole - <B> not found."));
    else
    {
        double x, y, z;
        x = parser->parseDouble(field.attribute("x"));
        y = parser->parseDouble(field.attribute("y"));
        z = parser->parseDouble(field.attribute("z"));
        B_val = Vector(x,y,z);
    }
    pugi::xml_node pos1 = node.child("p1");
    if (!pos1)
        throw(IOexception("InputParser::SoftEdgeDipole - <p1> not found."));
    else
    {
        double x, y, z;
        x = parser->parseDouble(pos1.attribute("x"));
        y = parser->parseDouble(pos1.attribute("y"));
        z = parser->parseDouble(pos1.attribute("z"));
        p_1 = Vector(x,y,z);
    }
    pugi::xml_node norm1 = node.child("n1");
    if (!norm1)
        throw(IOexception("InputParser::SoftEdgeDipole - <n1> not found."));
    else
    {
        double x, y, z;
        x = parser->parseDouble(norm1.attribute("x"));
        y = parser->parseDouble(norm1.attribute("y"));
        z = parser->parseDouble(norm1.attribute("z"));
        n_1 = Vector(x,y,z);
    }
    pugi::xml_node pos2 = node.child("p2");
    if (!pos2)
        throw(IOexception("InputParser::SoftEdgeDipole - <p2> not found."));
    else
    {
        double x, y, z;
        x = parser->parseDouble(pos2.attribute("x"));
        y = parser->parseDouble(pos2.attribute("y"));
        z = parser->parseDouble(pos2.attribute("z"));
        p_2 = Vector(x,y,z);
    }
    pugi::xml_node norm2 = node.child("n2");
    if (!norm2)
        throw(IOexception("InputParser::SoftEdgeDipole - <n2> not found."));
    else
    {
        double x, y, z;
        x = parser->parseDouble(norm2.attribute("x"));
        y = parser->parseDouble(norm2.attribute("y"));
        z = parser->parseDouble(norm2.attribute("z"));
        n_2 = Vector(x,y,z);
    }
    pugi::xml_node trans = node.child("transition");
    if (!trans)
        throw(IOexception("InputParser::SoftEdgeDipole - <transition> not found."));
    else
    {
        t_1 = parser->parseDouble(trans.attribute("l1"));
        t_2 = parser->parseDouble(trans.attribute("l2"));
    }
}

ElMagField SoftEdgeDipole::LocalField(double t, Vector X)
{
    Vector E = Vector(0.0, 0.0, 0.0);
    Vector B = Vector(0.0, 0.0, 0.0);
    // check distances from the limiting planes (outside positive)
    double d1 = dot(X-p_1,n_1);
    double d2 = dot(X-p_2,n_2);
    // scale the field at the edges
    double at1 = atan(2.0*d1/t_1);
    double at2 = atan(2.0*d2/t_2);
    double scal = 1.0 - (at1/Pi+0.5) - (at2/Pi+0.5);
    B = B_val * scal;
    return ElMagField(E, B);
}

