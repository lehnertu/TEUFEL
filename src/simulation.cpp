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
}

Simulation::~Simulation()
{
}

int Simulation::parseLattice(Lattice *lattice)
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
            type = element.attribute("type").value();
            name = element.attribute("name").value();
            if (type == "planar")
            {
                std::cout << name << "::PlanarUndulator" << std::endl;
                PlanarUndulator* Undu = new PlanarUndulator(element);
                lattice->addElement(Undu);
                /*
                std::cout << std::fixed << std::setprecision(3);
                std::cout << "  x=" << x << " m ";
                std::cout << "  y=" << y << " m ";
                std::cout << "  z=" << z << " m " << std::endl;
                std::cout << std::fixed << std::setprecision(6);
                std::cout << "  B= " << B << " T ";
                std::cout << "  lambda= " << period << " m ";
                std::cout << "  N= " << N << std::endl;
                */
            }
            else throw std::invalid_argument("unknown undulator type");
        }
        else throw std::invalid_argument("unknown lattice element");
    }
    return count;
}

void Simulation::parseBeam()
{
    pugi::xml_node beam = root.child("beam");
    if (!beam)
    {
	std::cout << "fatal error : node <beam> not found" << std::endl;
	exit(-1);
    };
}

void Simulation::parseObservations()
{
}
