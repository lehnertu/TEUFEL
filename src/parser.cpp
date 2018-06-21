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
#include <math.h>

#include "fields.h"
#include "global.h"
#include "parser.h"
#include "undulator.h"
#include "vector.h"

InputParser::InputParser(const pugi::xml_node node)
{
    // this is the <teufel> node
    root = node;
    calc = new mu::Parser();
    // we pre-define a number of constants in the calculator
    calc->DefineConst("_c", (double)SpeedOfLight);
    calc->DefineConst("_e", (double)ElementaryCharge);
    calc->DefineConst("_mec2", (double)mecsquared);
    calc->DefineConst("_eps0", (double)EpsNull);
    calc->DefineConst("_pi", (double)Pi);
    // evaluate all <calc> nodes present as children of root
    parseCalcChildren(root);
}

InputParser::~InputParser()
{
    delete calc;
}

void InputParser::parseCalc(const pugi::xml_node node)
{
    std::string prt = node.attribute("print").value();
    std::string var = node.attribute("var").value();
    std::string eq = node.attribute("eq").value();
    // if this is a print statement
    if (prt.length() > 0)
    {
        try
        {
            calc->SetExpr(eq);
            std::cout << "calc : " << prt << calc->Eval() << std::endl;
        }
        catch (mu::Parser::exception_type &e)
        {
            std::cout << "InputParser::parseCalc : can't evaluate " << eq << std::endl;
            std::cout << e.GetMsg() << endl;
        }
    };
    // if this is a variable definition
    if (var.length() > 0)
    {
        try
        {
            calc->SetExpr(eq);
            calc->DefineConst(var, calc->Eval());
        }
        catch (mu::Parser::exception_type &e)
        {
            std::cout << "InputParser::parseCalc : can't evaluate " << eq << std::endl;
            std::cout << e.GetMsg() << endl;
        }
    };
}

void InputParser::parseCalcChildren(const pugi::xml_node node)
{
    for (pugi::xml_node_iterator it = node.begin(); it != node.end(); ++it)
    {
        pugi::xml_node element = *it;
        std::string type = element.name();
        if ( type == "calc") parseCalc(element);
    }
}

double InputParser::parseValue(const pugi::xml_attribute attr)
{
    double retval = 0.0;
    try
    {
        calc->SetExpr(attr.value());
        retval = calc->Eval();
    }
    catch (mu::Parser::exception_type &e)
    {
        std::cout << "InputParser::parseValue : can't evaluate " << attr.value() << std::endl;
        std::cout << e.GetMsg() << endl;
    }
    return retval;
}

int InputParser::parseLattice(Lattice *lattice)
{
    int count = 0;
    std::string type;
    std::string name;
    pugi::xml_node latticenode = root.child("lattice");
    if (!latticenode)
        throw(IOexception("InputParser::parseLattice(Lattice - section <lattice> not found in input."));
    // loop over all children of the lattice node
    for (pugi::xml_node_iterator it = latticenode.begin(); it != latticenode.end(); ++it)
    {
        pugi::xml_node element = *it;
        type = element.name();
        if (type == "calc")
            parseCalc(element);
        else if (type == "undulator")
        {
            type = element.attribute("type").value();
            name = element.attribute("name").value();
            if (type == "planar")
            {
                std::cout << name << "::PlanarUndulator" << std::endl;
                // the undulator object parses its own input
                // we provide a reference to the parser
                PlanarUndulator* Undu = new PlanarUndulator(element, this);
                lattice->addElement(Undu);
            }
            else
                throw(IOexception("InputParser::parseLattice(Lattice - unknown undulator type."));
            count++;
        }
        else
            throw(IOexception("InputParser::parseLattice(Lattice - unknown lattice element."));
    }
    return count;
}

int InputParser::parseBeam(Beam *beam)
{
    int count = 0;
    pugi::xml_node beamnode = root.child("beam");
    if (!beamnode)
    {
        throw(IOexception("InputParser::parseBeam - section <beam> not found in input."));
    }
    else
    {
        for (pugi::xml_node_iterator it = beamnode.begin(); it != beamnode.end(); ++it)
        {
            pugi::xml_node entry = *it;
            std::string type = entry.name();
            if (type == "calc")
                parseCalc(entry);
            else if (type == "particle")
            {
                double gamma = parseValue(entry.attribute("gamma"));
                double betagamma = sqrt(gamma * gamma - 1.0);
                double charge = parseValue(entry.attribute("charge"))/ElementaryCharge;
                double cmr = parseValue(entry.attribute("cmr"));
                parseCalcChildren(entry);
                pugi::xml_node posnode = entry.child("position");
                if (!posnode)
                    throw(IOexception("InputParser::parseBeam - <particle> <position> not found."));
                double x, y, z;
                x = parseValue(posnode.attribute("x"));
                y = parseValue(posnode.attribute("y"));
                z = parseValue(posnode.attribute("z"));
                Vector pos = Vector(x, y, z);
                pugi::xml_node dirnode = entry.child("direction");
                if (!dirnode)
                    throw(IOexception("InputParser::parseBeam - <particle> <direction> not found."));
                x = parseValue(dirnode.attribute("x"));
                y = parseValue(dirnode.attribute("y"));
                z = parseValue(dirnode.attribute("z"));
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
            else
                throw(IOexception("InputParser::parseBeam - unknown type of beam entry."));
        };
    };
    return count;
}

void InputParser::parseTracking(Beam *beam)
{
    pugi::xml_node track = root.child("tracking");
    if (!track)
    {
        throw(IOexception("InputParser::parseTracking - section <tracking> not found in input."));
    }
    else
    {
        std::string method = track.attribute("method").value();
        if (method == "Vay")
        {
            beam->setTrackingMethod(TRACKING_VAY);
        }
        else
            throw(IOexception("InputParser::parseTracking - unknown tracking method."));
        pugi::xml_attribute timestep = track.attribute("delta_t");
        if (!timestep)
            throw(IOexception("InputParser::parseTracking - <tracking> attribute delta_t not found."));
        beam->setTimeStep(parseValue(timestep));
        pugi::xml_attribute nst = track.attribute("n");
        if (!nst)
            throw(IOexception("InputParser::parseTracking - <tracking> attribute n not found."));
        beam->setNOTS(nst.as_int());
    }
}

int InputParser::parseObservers(std::vector<Observer*> *listObservers, Beam *beam)
{
    int NOO = 0;
    pugi::xml_node obsroot = root.child("observer");
    if (!obsroot)
    {
        throw(IOexception("InputParser::parseObservers - section <observer> not found in input."));
    }
    else
    {
        for (pugi::xml_node_iterator it = obsroot.begin(); it != obsroot.end(); ++it)
        {
            pugi::xml_node obs = *it;
            std::string type = obs.name();
            if (type == "calc")
                parseCalc(obs);
            else if (type == "point")
            {
                NOO++;
            }
            else if (type == "screen")
            {
                double x, y, z;
                double t0, dt;
                pugi::xml_attribute fn = obs.attribute("file");
                if (!fn)
                    throw(IOexception("InputParser::parseObservers - <screen> attribute file not found."));                
                parseCalcChildren(obs);
                pugi::xml_node posnode = obs.child("position");
                if (!posnode)
                    throw(IOexception("InputParser::parseObservers - <screen> <position> not found."));
                x = parseValue(posnode.attribute("x"));
                y = parseValue(posnode.attribute("y"));
                z = parseValue(posnode.attribute("z"));
                Vector pos = Vector(x, y, z);
                pugi::xml_node hnode = obs.child("hpitch");
                if (!hnode)
                    throw(IOexception("InputParser::parseObservers - <screen> <hpitch> not found."));
                x = parseValue(hnode.attribute("x"));
                y = parseValue(hnode.attribute("y"));
                z = parseValue(hnode.attribute("z"));
                Vector dx = Vector(x, y, z);
                pugi::xml_attribute nh = hnode.attribute("n");
                if (!nh)
                    throw(IOexception("InputParser::parseObservers - <screen> <hpitch> attribute n not found."));
                pugi::xml_node vnode = obs.child("vpitch");
                if (!vnode)
                    throw(IOexception("InputParser::parseObservers - <screen> <vpitch> not found."));
                x = parseValue(vnode.attribute("x"));
                y = parseValue(vnode.attribute("y"));
                z = parseValue(vnode.attribute("z"));
                Vector dy = Vector(x, y, z);
                pugi::xml_attribute nv = vnode.attribute("n");
                if (!nv)
                    throw(IOexception("InputParser::parseObservers - <screen> <vpitch> attribute n not found."));
                pugi::xml_node tnode = obs.child("time");
                if (!tnode)
                    throw(IOexception("InputParser::parseObservers - <screen> <time> not found."));
                t0 = parseValue(tnode.attribute("t0"));
                dt = parseValue(tnode.attribute("dt"));
                pugi::xml_attribute nt = tnode.attribute("n");
                if (!nt)
                    throw(IOexception("InputParser::parseObservers - <screen> <time> attribute n not found."));
                ScreenObserver<Beam> *screenObs = new ScreenObserver<Beam>(
                    beam, fn.as_string(),
                    pos, dx, dy,
                    nh.as_int(),
                    nv.as_int(),
                    t0, dt, nt.as_int() );
                listObservers->push_back(screenObs);
                NOO++;
            }
            else
                throw(IOexception("InputParser::parseObservers - unknown observer type."));
        };
    }
    return NOO;
}
