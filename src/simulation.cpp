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

#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "fields.h"
#include "global.h"
#include "simulation.h"
#include "undulator.h"
#include "vector.h"

Simulation::Simulation(const pugi::xml_node node)
{
    root = node;
    lattice = new Lattice;
}

Simulation::~Simulation()
{
    delete lattice;
}

int Simulation::parseLattice()
{
    int count = 0;
    std::string type;
    std::string name;
    pugi::xml_node latticenode = root.child("lattice");
    if (!latticenode) throw std::invalid_argument("<lattice> not found");
    for (pugi::xml_node_iterator it = latticenode.begin(); it != latticenode.end(); ++it)
    {
	count++;
	pugi::xml_node element = *it;
	type = element.name();
	std::cout << "lattice::" << type << std::endl;
	if (type == "undulator")
	{
	    double x, y, z, B, period;
	    int N;
	    type = element.attribute("type").value();
	    name = element.attribute("name").value();
	    if (type == "planar")
	    {
		std::cout << name << "::PlanarUndulator" << std::endl;
		pugi::xml_node position = element.child("position");
		if (!position) throw std::invalid_argument("undulator <position> not found");
		else
		{
		    x = position.attribute("x").as_double(0.0);
		    y = position.attribute("y").as_double(0.0);
		    z = position.attribute("z").as_double(0.0);
		}
		pugi::xml_node field = element.child("field");
		if (!field) throw std::invalid_argument("undulator <field> not found");
		else
		{
		    B = field.attribute("B").as_double(0.0);
		    period = field.attribute("period").as_double(0.0);
		    N = field.attribute("N").as_int(0);
		}
		std::cout << std::fixed << std::setprecision(3);
		std::cout << "  x=" << x << " m ";
		std::cout << "  y=" << y << " m ";
		std::cout << "  z=" << z << " m " << std::endl;
		std::cout << std::fixed << std::setprecision(6);
		std::cout << "  B= " << B << " T ";
		std::cout << "  lambda= " << period << " m ";
		std::cout << "  N= " << N << std::endl;
		// now we can construct the undulator
		PlanarUndulator* Undu = new PlanarUndulator(Vector(x,y,z));
		Undu->Setup(B, period, N);
		lattice->addElement(Undu);
	    } else throw std::invalid_argument("unknown undulator type");
	}
	else throw std::invalid_argument("unknown lattice element");
		std::string type = element.name();
	std::cout << std::endl;
    }
    return count;
}

int Simulation::parseBeam()
{
    int count = 0;
    pugi::xml_node beam = root.child("beam");
    if (!beam)
    {
	std::cout << "fatal error : node <beam> not found" << std::endl;
	exit(-1);
    };
    return count;
}

void Simulation::run()
{
}

void Simulation::generateOutput()
{
}
