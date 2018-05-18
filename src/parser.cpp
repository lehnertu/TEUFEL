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
#include <math.h>

#include "fields.h"
#include "global.h"
#include "parser.h"
#include "undulator.h"
#include "vector.h"

InputParser::InputParser(const pugi::xml_node node)
{
    root = node;
}

InputParser::~InputParser()
{
}

int InputParser::parseLattice(Lattice *lattice)
{
    int count = 0;
    std::string type;
    std::string name;
    pugi::xml_node latticenode = root.child("lattice");
    if (!latticenode) throw std::invalid_argument("section <lattice> not found in input.");
    // loop over all children of the lattice node
    for (pugi::xml_node_iterator it = latticenode.begin(); it != latticenode.end(); ++it)
    {
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
                // the undulator object parses its own input
                PlanarUndulator* Undu = new PlanarUndulator(element);
                lattice->addElement(Undu);
            }
            else throw std::invalid_argument("unknown undulator type");
            count++;
        }
        else throw std::invalid_argument("unknown lattice element");
    }
    return count;
}

int InputParser::parseBeam(Beam *beam)
{
    int count = 0;
    pugi::xml_node beamnode = root.child("beam");
    if (!beamnode)
    {
        throw std::invalid_argument("section <beam> not found in input.");
    }
    else
    {
        // count the number of objects found within the beam
        for (pugi::xml_node_iterator it = beamnode.begin(); it != beamnode.end(); ++it)
        {
            pugi::xml_node entry = *it;
            std::string type = entry.name();
            if (type == "particle")
            {
                double gamma = entry.attribute("gamma").as_double(0.0);
                double betagamma = sqrt(gamma * gamma - 1.0);
                double charge = entry.attribute("charge").as_double(0.0)/ElementaryCharge;
                double cmr = entry.attribute("cmr").as_double(0.0);
                pugi::xml_node posnode = entry.child("position");
                if (!posnode) throw std::invalid_argument("<particle> <position> not found");
                double x, y, z;
                x = posnode.attribute("x").as_double(0.0);
                y = posnode.attribute("y").as_double(0.0);
                z = posnode.attribute("z").as_double(0.0);
                Vector pos = Vector(x, y, z);
                pugi::xml_node dirnode = entry.child("direction");
                if (!dirnode) throw std::invalid_argument("<particle> <direction> not found");
                x = dirnode.attribute("x").as_double(0.0);
                y = dirnode.attribute("y").as_double(0.0);
                z = dirnode.attribute("z").as_double(0.0);
                Vector mom = Vector(x, y, z);
                mom.normalize();
                mom *= betagamma;
                Vector acc = Vector(0.0, 0.0, 0.0);
                // now we have all information - create a particle
                ChargedParticle *p = new ChargedParticle(charge,cmr*charge);
                p->initTrajectory(0.0, pos, mom, acc);
                // create a bunch with just this particle and add it to the beam
                Bunch *single = new Bunch();
                single->Add(p);
                beam->Add(single);
                count++;
            }
            else throw std::invalid_argument("unknown type of beam entry");
        };
    };
    return count;
}

void InputParser::parseObservations()
{
}
