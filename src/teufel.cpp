/*=========================================================================

  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers

  Module:    undulator test case

  Copyright (c) 2018 U. Lehnert

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

/*!
    \brief Radiation from the ELBE U300 THz source

    @author Ulf Lehnert
    @date 16.3.2018
    @file teufel.cpp
    
    This test case tracks a number of particle in an undulator field.
    The particle energy is 24 MeV. we are tracking a number of macroparticles
    distributed as a number of bunches over the compute nodes. The total beam
    corresponds to 4.37e8 electrons that is a charge of 70 pC.
    As there is no interaction between the particles all bunches are
    tracked independently in parallel.
    
    The particles propagates through an undulator of 8 periods with 300 mm
    preiod length. The particles start at z=0, the undulator is centered at z=2.0m.
    At z=10m the produced radiation is observed.
    
    At c*t=2.0m a snapshot of the particle distribution and the local fields is generated.

*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "config.h"
#include "fields.h"
#include "global.h"
#include "observer.h"
#include "parser.h"

#include "pugixml.hpp"

int main(int argc, char* argv[])
{
    std::cout << std::endl;
    std::cout << "TEUFEL " << TEUFEL_VERSION_MAJOR << ".";
    std::cout.width(2);
    std::cout.fill('0');
    std::cout << TEUFEL_VERSION_MINOR << ".";
    std::cout.width(2);
    std::cout.fill('0');
    std::cout << TEUFEL_VERSION_PATCH << std::endl;
    cout << std::endl <<"THz-Emission From Undulators and Free-Electron Lasers" << std::endl << std::endl;
    
    // The first command line argument is interpreted as the input file name
    // This is an XML document which is opened and parsed here.
    pugi::xml_document doc;
    if (argc < 2) {
        std::cout << "Usage is teufel <infile>" << std::endl << std::endl;
        exit(1);
    } else {
        std::cout << "reading XML input from " << argv[1] << std::endl;
        pugi::xml_parse_result result = doc.load_file(argv[1]);
        if (result)
        {
            std::cout << "input parsed without errors" << std::endl;
        }
        else
        {
            std::cout << "ERROR reading file " << argv[1] << std::endl;
            std::cout << "Error description: " << result.description() << std::endl;
            exit(1);
        }
    };
    pugi::xml_node root = doc.child("teufel");
    if (!root)
        throw(IOexception("TEUFEL::InputParser - root node <teufel> not found."));
    string description = root.attribute("description").value();
    string author = root.attribute("author").value();
    std::cout << "case : " << description << std::endl;
    std::cout << "by : " << author << std::endl << std::endl;
    
    // Further parsing of the input document is done by the Parser object
    InputParser *parse = new InputParser(root);
    
    // We create an empty lattice object.
    // All lattice elements found when parsing the input are added to this.
    Lattice *lattice = new Lattice;
    int NoE = parse->parseLattice(lattice);
    std::cout << std::endl;
    std::cout << "lattice of " << NoE << " elements created." << std::endl;
    
    // We create an empty beam object.
    // Then we call the parser to fill in the necessary information
    // from the input file.
    Beam *beam = new Beam();
    int NoB = parse->parseBeam(beam);
    std::cout << std::endl;
    std::cout << "beam of " << NoB << " bunches created." << std::endl;
    std::cout << "total number of particles : " << beam->getNOP() << std::endl;
    std::cout << "total charge : " << beam->getTotalCharge()*ElementaryCharge*1.0e9 << "nC" << std::endl;
    std::cout << std::endl;

    // get all tracking information from the input file
    std::vector<watch_t> watches;
    parse->parseTracking(beam, &watches);
    if (teufel::rank==0) std::cout << "defined " << (int)watches.size() << " watch points." << std::endl;
    if (teufel::rank==0) std::cout << std::endl;

    // parse all observer definitions
    std::vector<Observer*> listObservers;
    parse->parseObservers(&listObservers);
    std::cout << std::endl;
    std::cout << "defined " << (int)listObservers.size() << " observers." << std::endl;
    std::cout << std::endl;
    
    // we are done with the input document
    delete parse;
    
    // prepare the tracking of the beam
    beam->setupTracking(lattice);

    // handle watch point of initial particle distribution
    for (int iw=0; iw<(int)watches.size(); iw++)
    {
        watch_t w = watches.at(iw);
        if (w.step == 0)
        {
            std::cout << "writing watch point " << w.filename << std::endl;
            int nw = beam->WriteWatchPointHDF5(w.filename);
            std::cout << nw << " particles written." << std::endl;
        }
    }
    // do the tracking of the beam
    std::cout << std::endl;
    std::cout << "tracking " << beam->getNOP() << " particles ..." << std::endl;
    // record the start time
    timespec start_time, current_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);
    timespec print_time = start_time;
    for (int step=0; step<beam->getNOTS(); step++)
    {
        beam->doStep(lattice);
        // handle watch points
        for (int iw=0; iw<(int)watches.size(); iw++)
        {
            watch_t w = watches.at(iw);
            if (w.step == step+1)
            {
                std::cout << "writing watch point " << w.filename << std::endl;
                int nw = beam->WriteWatchPointHDF5(w.filename);
                if (nw>=0)
                {
                    std::cout << nw << " particles written." << std::endl;
                    std::cout << std::endl;
                }
                else
                    std::cout << "error writing " << w.filename << std::endl;
            }
        }
        // make a print once every 10s
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &current_time);
        double elapsed = current_time.tv_sec-print_time.tv_sec +
            1e-9*(current_time.tv_nsec-print_time.tv_nsec);
        if (elapsed>10.0)
        {
            print_time = current_time;
            std::cout << "tracking step " << step << std::endl;
        };
    }
    // record the finish time
    timespec stop_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop_time);
    std::cout << "finished tracking particles." << std::endl;
    double elapsed = stop_time.tv_sec-start_time.tv_sec +
        1e-9*(stop_time.tv_nsec-start_time.tv_nsec);
    std::cout << "time elapsed during tracking : " << elapsed << " s" << std::endl;

    // compute all observations
    for (int i=0; i<(int)listObservers.size(); i++)
    {
        std::cout << std::endl << "computing observer No. " << i+1 << std::endl;
        Observer *obs = listObservers.at(i);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);
        std::cout << "integrating ... " << std::endl;
        obs->integrate(beam);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop_time);
        double elapsed = stop_time.tv_sec-start_time.tv_sec +
            1e-9*(stop_time.tv_nsec-start_time.tv_nsec);
        std::cout << "time elapsed : " << elapsed << " s" << std::endl;
        obs->generateOutput();
    }
        
    // delete all observers
    for (int i=0; i<(int)listObservers.size(); i++)
        delete listObservers.at(i);
    listObservers.clear();
    
    // deleting the beam in recursion deletes all contained bunches and particles
    delete beam;
    
    // deleting the lattice automatically deletes all lattice elments
    delete lattice;
    
    return 0;
}
