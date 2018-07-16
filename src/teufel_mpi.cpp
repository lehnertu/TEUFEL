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
    \brief TEUFEL main program

    @author Ulf Lehnert
    @date 12.7.2018
    @file teufel_mpi.cpp
    
    This is the MPI version of main executable of TEUFEL.
    It should be executed like
    
    mpiexec -n 4 teufel_mpi input.xml
    
    All computation is distributed over a number of nodes connected
    in a CPU cluster. Even on sigle-processor machines it makes sense to
    use the MPI parallelization to equally distribute the load over
    all cores available.
    
    The input file is read and interpreted by all nodes individually.
    The beam is only created on the root node and ditributed over all
    nodes such that all receive approximately equal numbers of particles.
    Then the tracking is done by all nodes independently for their own
    set of particles.
    
*/

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <time.h>

#include "bunch.h"
#include "beam.h"
#include "config.h"
#include "fields.h"
#include "global.h"
#include "observer.h"
#include "parser.h"
#include "particle.h"
#include "undulator.h"

#include "pugixml.hpp"

int NOP = 1e2;          // number of particles
int NOTS = 4000;        // number of time steps

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);
    
    int NumberOfCores = 1;
    teufel::rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &NumberOfCores);
    MPI_Comm_rank(MPI_COMM_WORLD, &teufel::rank);
    
    if (teufel::rank==0)
    {
        std::cout << std::endl;
        std::cout << " TEUFEL " << TEUFEL_VERSION_MAJOR << ".";
        std::cout.width(2);
        std::cout.fill('0');
        std::cout << TEUFEL_VERSION_MINOR << ".";
        std::cout.width(2);
        std::cout.fill('0');
        std::cout << TEUFEL_VERSION_PATCH << std::endl;
        cout << std::endl <<" THz-Emission From Undulators and Free-Electron Lasers" << std::endl << std::endl;
        cout << std::endl <<" TEUFEL parallel computing on " << NumberOfCores << " cores." << std::endl << std::endl;
    }

    // ===============================================================
    // parsing the input file by all processes in parallel
    // console output is only printed from rank 0
    // ===============================================================
    
    // The first command line argument is interpreted as the input file name.
    // This is an XML document which is opened and parsed here.
    pugi::xml_document doc;
    if (argc < 2) {
        if (teufel::rank==0) std::cout << "Usage is teufel <infile>" << std::endl << std::endl;
        exit(1);
    } else {
        if (teufel::rank==0) std::cout << "reading XML input from " << argv[1] << std::endl;
        pugi::xml_parse_result result = doc.load_file(argv[1]);
        if (result)
        {
            if (teufel::rank==0) std::cout << "input parsed without errors" << std::endl;
        }
        else
        {
            if (teufel::rank==0) std::cout << "ERROR reading file " << argv[1] << std::endl;
            if (teufel::rank==0) std::cout << "Error description: " << result.description() << std::endl;
            exit(1);
        }
    };
    pugi::xml_node root = doc.child("teufel");
    if (!root)
        throw(IOexception("TEUFEL::InputParser - root node <teufel> not found."));
    string description = root.attribute("description").value();
    string author = root.attribute("author").value();
    if (teufel::rank==0) std::cout << "case : " << description << std::endl;
    if (teufel::rank==0) std::cout << "by : " << author << std::endl << std::endl;
    
    // Further parsing of the input document is done by the parser object.
    InputParser *parse = new InputParser(root);
    
    // We create an empty lattice object.
    // All lattice elements found when parsing the input are added to this.
    Lattice *lattice = new Lattice;
    int NoE = parse->parseLattice(lattice);
    if (teufel::rank==0) std::cout << std::endl;
    if (teufel::rank==0) std::cout << "lattice of " << NoE << " elements created." << std::endl;
    
    // We create an empty beam object.
    // Then we call the parser to fill in the necessary information
    // from the input file.
    Beam *beam = new Beam();
    int NoB = parse->parseBeam(beam);
    if (teufel::rank==0) std::cout << std::endl;
    if (teufel::rank==0) std::cout << "beam of " << NoB << " bunches created." << std::endl;
    if (teufel::rank==0) std::cout << "total number of particles : " << beam->getNOP() << std::endl;
    if (teufel::rank==0) std::cout << "total charge : " << beam->getTotalCharge()*ElementaryCharge*1.0e9 << "nC" << std::endl;
    if (teufel::rank==0) std::cout << std::endl;

    // get all tracking information from the input file
    std::vector<watch_t> watches;
    int NoW = parse->parseTracking(beam, &watches);
    if (teufel::rank==0) std::cout << "defined " << NoW << " watch points." << std::endl;
    if (NoW != (int)watches.size())
        if (teufel::rank==0) std::cout << "WARNING : number of watches (" << NoW << ") differs from length of the list (" << watches.size() << ")" << std::endl;
    if (teufel::rank==0) std::cout << std::endl;
    
    // parse all observer definitions
    std::vector<Observer*> listObservers;
    int NoO = parse->parseObservers(&listObservers);
    if (teufel::rank==0) std::cout << std::endl;
    if (teufel::rank==0) std::cout << "defined " << NoO << " observers." << std::endl;
    if (NoO != (int)listObservers.size())
        if (teufel::rank==0) std::cout << "WARNING : number of observers (" << NoO << ") differs from length of the list (" << listObservers.size() << ")" << std::endl;
    if (teufel::rank==0) std::cout << std::endl;
    
    // done parsing the input file
    delete parse;

    // ===============================================================
    // re-distribute all particles into one bunch per node
    // all containing approximately equal numbers of particles
    // ===============================================================

    // buffers for MPI communication
    double *sendbuffer = new double[PARTICLE_SERIALIZE_BUFSIZE];
    double *recbuffer = new double[PARTICLE_SERIALIZE_BUFSIZE];
    MPI_Request send_req;

    // we keep track of the particle numbers per node
    int pc_total = 0;
    int *pc_node = new int[NumberOfCores];
    
    Bunch *trackedBunch = new Bunch();
    if (teufel::rank == 0) std::cout << "distribute particles..." << std::endl;
    // all nodes synchronously traverse the beam
    for (int i=0; i<NumberOfCores; i++) pc_node[i]=0;
    int nob = beam->getNOB();
    for (int i=0; i<nob; i++)
    {
        Bunch *b = beam->getBunch(i);
        int nop = b->getNOP();
        for (int j=0; j<nop; j++)
        {
            // only the root node is sending particle information
            if (teufel::rank == 0)
            {
                ChargedParticle *p = b->getParticle(j);
                p->serialize(sendbuffer);
                // non-blocking send
                MPI_Issend(sendbuffer, PARTICLE_SERIALIZE_BUFSIZE, MPI_DOUBLE,
                    pc_total % NumberOfCores, 42, MPI_COMM_WORLD, &send_req);
            };
            // all nodes in turn receive the particles
            if (teufel::rank == (pc_total % NumberOfCores) )
            {
                MPI_Recv(recbuffer, PARTICLE_SERIALIZE_BUFSIZE, MPI_DOUBLE,
                    0, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                ChargedParticle *tp = new ChargedParticle(recbuffer);
                trackedBunch->Add(tp);
            };
            // the root node waits for completion of the transfer
            if (teufel::rank == 0)
            {
                MPI_Wait(&send_req,MPI_STATUS_IGNORE);
            };
            pc_node[pc_total % NumberOfCores]++;
            pc_total++;
        };
    };
    
    MPI_Barrier(MPI_COMM_WORLD);

    // ===============================================================
    // track the beam on all nodes in parallel
    // ===============================================================

    // this is the beam we will actually track
    Beam *trackedBeam = new Beam();
    trackedBeam->Add(trackedBunch);
    // copy the tracking information
    trackedBeam->setTrackingMethod(beam->getTrackingMethod());
    trackedBeam->setTimeStep(beam->getTimeStep());
    trackedBeam->setNOTS(beam->getNOTS());
    // prepare the tracking of the beam
    trackedBeam->setupTracking(lattice);
    
    if (trackedBeam->getNOP() != pc_node[teufel::rank])
    {
        std::cout << "node " << teufel::rank << " has wrong number of particles - should be" << pc_node[teufel::rank] << std::endl;
        throw(IOexception("TEUFEL internal error - aborting."));
    };
    std::cout << "node " << teufel::rank << " tracking beam of " <<
        trackedBeam->getNOP() << " particles." << std::endl;
    
    MPI_Barrier(MPI_COMM_WORLD);

    // handle watch point of initial particle distribution if requested
    for (int iw=0; iw<NoW; iw++)
    {
        watch_t w = watches.at(iw);
        if (w.step == 0)
        {
            if (teufel::rank == 0)
                std::cout << "gathering data for watch point " << w.filename << std::endl;
            // we use the initial beam object to gather the particle information
            beam->clear();
            Bunch *outb = new Bunch();
            Bunch *b = trackedBeam->getBunch(0);
            // iterate over the nodes and particles
            for (int n=0; n<NumberOfCores; n++)
                for (int j=0; j<pc_node[n]; j++)
                {
                    if (teufel::rank == n)
                    {
                        ChargedParticle *p = b->getParticle(j);
                        p->serialize(sendbuffer);
                        MPI_Send(sendbuffer, PARTICLE_SERIALIZE_BUFSIZE, MPI_DOUBLE,
                                 0, 42, MPI_COMM_WORLD);
                    };
                    if (teufel::rank == 0)
                    {
                        MPI_Recv(recbuffer, PARTICLE_SERIALIZE_BUFSIZE, MPI_DOUBLE,
                                 n, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        ChargedParticle *outp = new ChargedParticle(recbuffer);
                        outb->Add(outp);
                    };
                };
            MPI_Barrier(MPI_COMM_WORLD);
            if (teufel::rank == 0)
            {
                std::cout << "writing watch point " << w.filename << std::endl;
                beam->Add(outb);
                int nw = beam->WriteWatchPointHDF5(w.filename);
                std::cout << nw << " particles written." << std::endl;
                std::cout << std::endl;
            };
            // this deletes all contained bunches and particles
            beam->clear();
        }
    }

    // record the start time
    double start_time = MPI_Wtime();
    double print_time = start_time;

    // do the tracking of the beam
    for (int step=0; step<trackedBeam->getNOTS(); step++)
    {
//         // do a step
        trackedBeam->doStep(lattice);
        // handle watch points
        for (int iw=0; iw<NoW; iw++)
        {
            watch_t w = watches.at(iw);
            if (w.step == step+1)
            {
                if (teufel::rank == 0)
                    std::cout << "gathering data for watch point " << w.filename << std::endl;
                // we use the initial beam object to gather the particle information
                beam->clear();
                Bunch *outb = new Bunch();
                Bunch *b = trackedBeam->getBunch(0);
                // iterate over the nodes and particles
                for (int n=0; n<NumberOfCores; n++)
                    for (int j=0; j<pc_node[n]; j++)
                    {
                        if (teufel::rank == n)
                        {
                            ChargedParticle *p = b->getParticle(j);
                            p->serialize(sendbuffer);
                            MPI_Send(sendbuffer, PARTICLE_SERIALIZE_BUFSIZE, MPI_DOUBLE,
                                    0, 42, MPI_COMM_WORLD);
                        };
                        if (teufel::rank == 0)
                        {
                            MPI_Recv(recbuffer, PARTICLE_SERIALIZE_BUFSIZE, MPI_DOUBLE,
                                    n, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            ChargedParticle *outp = new ChargedParticle(recbuffer);
                            outb->Add(outp);
                        };
                    };
                MPI_Barrier(MPI_COMM_WORLD);
                if (teufel::rank == 0)
                {
                    std::cout << "writing watch point " << w.filename << std::endl;
                    beam->Add(outb);
                    int nw = beam->WriteWatchPointHDF5(w.filename);
                    std::cout << nw << " particles written." << std::endl;
                    std::cout << std::endl;
                };
                // this deletes all contained bunches and particles
                beam->clear();
            }
        }

        // make a print once every 10s
        if (teufel::rank == 0)
        {
            double current_time = MPI_Wtime();
            if (current_time-print_time > 10)
            {
                print_time = current_time;
                std::cout << "tracking step " << step << std::endl;
            };
        };
    }

    // record the finish time
    if (teufel::rank == 0)
    {
        double stop_time = MPI_Wtime();
        std::cout << "finished tracking particles." << std::endl;
        std::cout << "time elapsed during tracking : " << stop_time-start_time << " s" << std::endl;
    };

    // ===============================================================
    // compute the radiation observations
    // each node does the computation for it's own set of particles
    // the fields are combined before writing the files
    // ===============================================================

    MPI_Barrier(MPI_COMM_WORLD);
    
    // compute all observations
    for (int i=0; i<NoO; i++)
    {
        double start_time = MPI_Wtime();
        if (teufel::rank == 0)
        {
            std::cout << std::endl << "computing observer No. " << i+1 << std::endl;
        };
        Observer *obs = listObservers.at(i);
        std::cout << "node " << teufel::rank << " integrating ... " << std::endl;
        obs->integrate(trackedBeam);
        double stop_time = MPI_Wtime();
        if (teufel::rank == 0)
        {
            std::cout << "time elapsed : " << stop_time-start_time << " s" << std::endl;
        };
        
        // collect all the field computed on the individual nodes into the master node
        unsigned int count = obs->getBufferSize();
        std::cout << "Node " << teufel::rank << " allocating buffers for "<< count << " doubles" << std::endl;
        // fill the buffer and get its address
        double* nodeBuffer = obs->getBuffer();
        if (nodeBuffer==0)
        {
            std::cout << "MPI_Reduce for Obsserver : node " << teufel::rank << " was unable to get the node buffer." << std::endl;
            throw(IOexception("memory allocation error"));
        };
        // we have to define a buffer for the sum on all nodes
        double* reduceBuffer = new double[count];
        if (reduceBuffer==0)
        {
            std::cout << "MPI_Reduce for Obsserver : node " << teufel::rank << " was unable to get the recude buffer." << std::endl;
            throw(IOexception("memory allocation error"));
        };
        MPI_Barrier(MPI_COMM_WORLD);
        if (teufel::rank == 0)
            std::cout << std::endl << "All buffers allocated." << std::endl;
        MPI_Reduce(
            nodeBuffer,                 // send buffer
            reduceBuffer,               // receive buffer
            count,                      // number of values
            MPI_DOUBLE,                 // data type
            MPI_SUM,                    // operation
            0,                          // rank of the root process
            MPI_COMM_WORLD              // communicator
        );
        if (teufel::rank == 0)
            std::cout << "Reduce finished." << std::endl;
        // the root node copies the data from the reduce buffer into the Observer
        // and writes it to the output file
        if (teufel::rank == 0)
        {
            obs->fromBuffer(reduceBuffer,count);
            try
            { 
                    obs->generateOutput();
                    std::cout << "Observer output file written." << std::endl;
            }
            catch (exception& e) { cout << e.what() << endl;}
        }
        // if (teufel::rank == 0) delete reduceBuffer;
        delete reduceBuffer;
        delete nodeBuffer;
    }
        

/*

    // collect all the field computed on the individual nodes
    // into the master node
*/

    // delete all observers
    for (int i=0; i<NoO; i++)
        delete listObservers.at(i);
    listObservers.clear();
    
    // delete all watches
    watches.clear();

    // deleting the lattice automatically deletes all lattice elments
    delete lattice;
    
    // deleting the beam automatically deletes all bunches and particles belonging to it
    delete beam;
    delete trackedBeam;

    delete sendbuffer;
    delete recbuffer;
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (teufel::rank==0)
        std::cout << std::endl << " TEUFEL run finished." << std::endl << std::endl;
    
    MPI_Finalize();
    return 0;
}
