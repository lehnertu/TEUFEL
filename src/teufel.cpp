/*=========================================================================

  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers

  Module:    Main Program

  Copyright (c) 2017 U. Lehnert

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/

#include "pugixml.hpp"
#include <iostream>

#include "global.h"
#include "simulation.h"

int main ( int argc, char *argv[] )
{
    std::cout << std::endl;
    std::cout << "TEUFEL - THz Emission from Undulators and Free-Electron Lasers" << std::endl;
    std::cout << "Ulf Lehnert 9/2017" << std::endl;
    std::cout << std::endl;

    // we should have exactly one argument - the input file name
    if ( argc != 2 )
    {
	cout<<"usage: "<< argv[0] <<" <filename>\n";
	exit(-1);
    };
    
    pugi::xml_document doc;
    if (!doc.load_file(argv[1]))
    {
	std::cout << "input file read error" << std::endl;
	exit(-1);
    };

    pugi::xml_node root = doc.child("teufel");
    if (!root)
    {
	std::cout << "fatal error : root node <teufel> not found" << std::endl;
	exit(-1);
    };

    std::cout << root.attribute("description").value() << std::endl;
    std::cout << root.attribute("author").value() << std::endl;
    std::cout << std::endl;
    
    Simulation *sim = new Simulation(root);
    int num = sim->parseLattice();
    std::cout << "parsed " << num << " lattice elements" << std::endl;
    num = sim->parseBeam();
    
    delete sim;
}

