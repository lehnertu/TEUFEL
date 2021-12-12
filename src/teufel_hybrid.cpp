/*=========================================================================

  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers

  Module:    main executable for particle tracking and radiation emission

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
    @date 14.7.2021
    @file teufel_hybrid.cpp
    
    This is the hybrid version of the main executable of TEUFEL.
    It should be executed like
    
    export OMP_NUM_THREADS=8
    mpiexec -n 4 teufel_hybrid input.xml
    
    All computation is distributed over a number of nodes connected
    in a CPU cluster. The communication between the nodes is done with MPI.
    On each node the load is distibuted over the available cores using OpenMP.
    
    The input file is read and interpreted by all nodes individually.
    The beam is only created on the root node and ditributed over all nodes.
    Each node holds the complete trajectory information of all particles
    but tracks only its own subset. After every tracking step the trajectory
    information is communicated between the nodes.
    
*/

#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
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

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);
    
    int NumberOfNodes = 1;
    teufel::rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &NumberOfNodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &teufel::rank);
    int NumberOfThreads = omp_get_max_threads();
    
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
        cout << std::endl <<" TEUFEL parallel computing on " << NumberOfNodes << " nodes." << std::endl;
        cout << std::endl;
    }
    cout << std::flush;
    usleep(100000);
    MPI_Barrier(MPI_COMM_WORLD);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int namelen;
    int mem_size = -1;
    int PID = getpid();
    MPI_Get_processor_name(processor_name, &namelen);
    std::string line;
    std::string delimiter = ":";
    std::ifstream myfile;
    std::stringstream ss;
    myfile.open("/proc/meminfo");
    if(!myfile.is_open()) {
        std::cout << "cannot access /proc/meminfo - memory size unknown." << std::endl;
    } else {
        while(getline(myfile, line)) {
            // std::cout << line << std::endl;
            size_t pos = line.find(delimiter);
            std::string name = line.substr(0, pos);
            std::string rest = line.erase(0,pos+1);
            if (name=="MemTotal") {
                mem_size = std::stoi(rest);
            }
        }
        myfile.close();
    }

    cout << "node " << teufel::rank << " : " << processor_name
        << " PID=" << PID
        << " total memory " << mem_size/1e3 << " MB"
        << " using " << NumberOfThreads << " parallel threads"
        << std::endl;

    ss.str(std::string());
    ss << "/proc/" << PID << "/status";
    myfile.open(ss.str());
    if(!myfile.is_open()) {
        cout << "cannot access /proc/{PID}/status - memory usage unknown." << std::endl;
    } else {
        while(getline(myfile, line)) {
            // std::cout << line << std::endl;
            size_t pos = line.find(delimiter);
            std::string name = line.substr(0, pos);
            std::string rest = line.erase(0,pos+1);
            if (name=="VmSize") {
                mem_size = std::stoi(rest);
            }
        }
        cout << "node " << teufel::rank << " memory usage : " << mem_size/1e3 << " MB" << std::endl;
        myfile.close();
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100000);
    if (teufel::rank==0) cout << std::endl << std::flush;
    usleep(100000);
    MPI_Barrier(MPI_COMM_WORLD);
    
    // ===============================================================
    // parsing the input file by all nodes in parallel
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
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100000);
    if (teufel::rank==0) cout << std::endl << std::flush;
    usleep(100000);
    MPI_Barrier(MPI_COMM_WORLD);
    
    // We create an empty beam object.
    // Then we call the parser to fill in the necessary information from the input file.
    // The structure of this beam shall then be preserved troughout the tracking procedure.
    // That is necessary to keep the references to the contained bunches for logging
    // or other bunch-related procedures. This way it is possible to log information
    // for inidividual bunches making up the total beam.
    // After every tracking step, the current particle information distributed over
    // the compute nodes will be re-gathered into this object on all nodes.
    Beam *masterBeam = new Beam();
    std::vector<TrackingLogger<Bunch>*> listLoggers;
    int NoB = parse->parseBeam(masterBeam, &listLoggers);
    if (teufel::rank==0) std::cout << std::endl;
    if (teufel::rank==0) std::cout << "beam of " << NoB << " bunches created." << std::endl;
    if (teufel::rank==0) std::cout << "total number of particles : " << masterBeam->getNOP() << std::endl;
    if (teufel::rank==0) std::cout << "total charge : "
        << masterBeam->getTotalCharge()*ElementaryCharge*1.0e9 << "nC" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100000);
    if (teufel::rank==0) cout << std::endl << std::flush;
    usleep(100000);
    MPI_Barrier(MPI_COMM_WORLD);

    // get all tracking information from the input file
    std::vector<watch_t> watches;
    parse->parseTracking(masterBeam, &watches);
    if (teufel::rank==0) std::cout << "defined " << (int)watches.size() << " watch points." << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100000);
    if (teufel::rank==0) std::cout << std::endl << std::flush;
    usleep(100000);
    MPI_Barrier(MPI_COMM_WORLD);
    
    // parse all observer definitions
    std::vector<Observer*> listObservers;
    parse->parseObservers(&listObservers);
    if (teufel::rank==0) std::cout << std::endl;
    if (teufel::rank==0) std::cout << "defined " << (int)listObservers.size() << " observers." << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100000);
    if (teufel::rank==0) std::cout << std::endl << std::flush;
    usleep(100000);
    MPI_Barrier(MPI_COMM_WORLD);

    // We register the number of time steps to track the beam.
    // This has been set by the parser, but will be cleared when
    // initalizing the trajectories for actual tracking.
    int NOTS = masterBeam->getNOTS();
    if (teufel::rank==0) std::cout << "tracking for " << NOTS << " time steps." << std::endl;
    if (teufel::rank==0) std::cout << "master beam trajectory storage : "
        << (double)masterBeam->getNOP() * NOTS * 10.0 * sizeof(double) / 1e6
        << " MB  for " << masterBeam->getNOP() << " particles" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100000);
    if (teufel::rank==0) std::cout << std::endl << std::flush;
    usleep(100000);
    MPI_Barrier(MPI_COMM_WORLD);

    // done parsing the input file
    delete parse;

    ss.str(std::string());
    ss << "/proc/" << PID << "/status";
    myfile.open(ss.str());
    if(!myfile.is_open()) {
        cout << "cannot access /proc/{PID}/status - memory usage unknown." << std::endl;
    } else {
        while(getline(myfile, line)) {
            // std::cout << line << std::endl;
            size_t pos = line.find(delimiter);
            std::string name = line.substr(0, pos);
            std::string rest = line.erase(0,pos+1);
            if (name=="VmSize") {
                mem_size = std::stoi(rest);
            }
        }
        cout << "node " << teufel::rank << " memory usage : " << mem_size/1e3 << " MB" << std::endl;
        myfile.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100000);
    if (teufel::rank==0) cout << std::endl << std::flush;
    usleep(100000);
    MPI_Barrier(MPI_COMM_WORLD);

    // ===============================================================
    // broadcast the masterBeam to all nodes
    // ===============================================================

    int MASTER_NOP = masterBeam->getNOP();
    int TRACKING_NOP = (MASTER_NOP + NumberOfNodes -1) / NumberOfNodes;
    if (teufel::rank==0) std::cout << "tracking " << TRACKING_NOP << " particles per node." << std::endl;
    if (teufel::rank==0) std::cout << "tracking beam trajectory storage : "
        << (double)TRACKING_NOP * NOTS * 10.0 * sizeof(double) / 1e6
        << " MB  for " << TRACKING_NOP << " particles" << std::endl;
    // allocate communication buffers holding the required number of particles
    // the master buffer must provide space for missing particles on some nodes
    int MASTER_BUFSIZE = TRACKING_NOP * NumberOfNodes * 10;
    int TRACKING_BUFSIZE = TRACKING_NOP * 10;
    if (MASTER_BUFSIZE < masterBeam->getStepBufferSize())
    {
        throw(IOexception("TEUFEL internal error: master buffer too small - aborting."));
        return(-1);
    };
    // buffer for MPI communication
    double *masterBuffer = new double[MASTER_BUFSIZE];
    double *trackingBuffer = new double[TRACKING_BUFSIZE];
    
    // fill the send buffer on the master node
    masterBeam->bufferStep(masterBuffer);

    // the buffer ins transfered from the master to all nodes
    MPI_Bcast(masterBuffer, MASTER_BUFSIZE, MPI_DOUBLE,
        0, MPI_COMM_WORLD);
        
    // update the beam from the buffer
    masterBeam->clearTrajectories();
    // we will use the master beam to store the gathered tracking data
    // so we pre-allocate the necessary trajectory storage for every particle
    masterBeam->preAllocate(NOTS+1);
    masterBeam->setStepFromBuffer(masterBuffer);
    // now all nodes have an identical beam containig all particles
    
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100000);
    if (teufel::rank==0) cout << "master beam setup done." << std::endl;
    if (teufel::rank==0) cout << std::endl << std::flush;
    usleep(100000);
    MPI_Barrier(MPI_COMM_WORLD);

    ss.str(std::string());
    ss << "/proc/" << PID << "/status";
    myfile.open(ss.str());
    if(!myfile.is_open()) {
        cout << "cannot access /proc/{PID}/status - memory usage unknown." << std::endl;
    } else {
        while(getline(myfile, line)) {
            // std::cout << line << std::endl;
            size_t pos = line.find(delimiter);
            std::string name = line.substr(0, pos);
            std::string rest = line.erase(0,pos+1);
            if (name=="VmSize") {
                mem_size = std::stoi(rest);
            }
        }
        cout << "node " << teufel::rank << " memory usage : " << mem_size/1e3 << " MB" << std::endl;
        myfile.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100000);
    if (teufel::rank==0) cout << std::endl << std::flush;
    usleep(100000);
    MPI_Barrier(MPI_COMM_WORLD);

    // ===============================================================
    // re-distribute all particles into one bunch per node
    // all containing approximately equal numbers of particles
    // (the last node may have less)
    // ===============================================================

    // this is the bunch every compute node will actually track
    Bunch *trackedBunch = new Bunch();
    if (teufel::rank == 0) std::cout << "distribute particles " << TRACKING_NOP << " per node" << std::endl;
    // particles are distributed in chunks of TRACKING_NOP to the individual nodes
    int current_node = 0;              // the node to which we assign particles
    int current_node_particles = 0;    // the number of particles already assigned to the current node
    // traverse the beam and copy particles
    int nob = masterBeam->getNOB();
    for (int i=0; i<nob; i++)
    {
        Bunch *b = masterBeam->getBunch(i);
        int nop = b->getNOP();
        for (int j=0; j<nop; j++)
        {
            ChargedParticle *master_particle = b->getParticle(j);
            // if the particle is for this node
            if (teufel::rank==current_node)
            {
                // create a new particle as an exact copy of the master beam particle
                ChargedParticle *p = new ChargedParticle(master_particle);
                p->preAllocate(NOTS+1);
                trackedBunch->Add(p);
            }
            current_node_particles++;
            // next node if we have enough particles here
            if (current_node_particles==TRACKING_NOP)
            {
                current_node_particles = 0;
                current_node++;
            }
        }
    }
    
    // ===============================================================
    // track the beam on all nodes in parallel
    // ===============================================================

    // this is the beam we will actually track (private per node)
    Beam *trackedBeam = new Beam();
    trackedBeam->Add(trackedBunch);
    // copy the tracking information from the master beam
    trackedBeam->setTrackingMethod(masterBeam->getTrackingMethod());
    trackedBeam->setTimeStep(masterBeam->getTimeStep());
    // prepare the tracking of the beam
    trackedBeam->setupTracking(lattice);
    
    std::cout << "node " << teufel::rank << " tracking a beam of " <<
        trackedBeam->getNOP() << " particles." << std::endl;
    ss.str(std::string());
    ss << "/proc/" << PID << "/status";
    myfile.open(ss.str());
    if(!myfile.is_open()) {
        cout << "cannot access /proc/{PID}/status - memory usage unknown." << std::endl;
    } else {
        while(getline(myfile, line)) {
            // std::cout << line << std::endl;
            size_t pos = line.find(delimiter);
            std::string name = line.substr(0, pos);
            std::string rest = line.erase(0,pos+1);
            if (name=="VmSize") {
                mem_size = std::stoi(rest);
            }
        }
        cout << "node " << teufel::rank << " memory usage : " << mem_size/1e3 << " MB" << std::endl;
        myfile.close();
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100000);
    if (teufel::rank==0) cout << std::endl << std::flush;
    usleep(100000);
    MPI_Barrier(MPI_COMM_WORLD);

    // handle watch point of initial particle distribution if requested
    if (teufel::rank == 0)
        for (int iw=0; iw<(int)watches.size(); iw++)
        {
            watch_t w = watches.at(iw);
            if (w.step == 0)
            {
                std::cout << "writing watch point " << w.filename << std::endl;
                int nw = masterBeam->WriteWatchPointHDF5(w.filename);
                std::cout << nw << " particles written." << std::endl;
                std::cout << std::endl;
            }
        }

    // update loggers with initial particle distribution
    if (teufel::rank == 0)
        for (int il=0; il<(int)listLoggers.size(); il++)
            listLoggers.at(il)->update();
                    
    // record the start time
    double start_time = MPI_Wtime();
    double print_time = start_time;
    double current_time = start_time;

    // do the tracking of the beam
    // Shared memory parallelization is used for field integration only.
    // It may be emploeyed in the future during tracking for
    // computing the interaction field between the particles.
    // With a-priori fields tracking will only be MPI parallelized.
    for (int step=0; step<NOTS; step++)
    {
    
        // do a step
        trackedBeam->doStep(lattice);

        // for (int i=0; i<trackedBunch->getNOP(); i++)
        // {
        //     ChargedParticle *p = trackedBunch->getParticle(i);
        //     p->StepVay(lattice);
        // }
            
        MPI_Barrier(MPI_COMM_WORLD);
        // distribute the step result to all nodes
        // each node buffers its own tracked bunch
        trackedBeam->bufferStep(trackingBuffer);
        // the buffers are gathered into the master buffer on all nodes
        MPI_Allgather(trackingBuffer, TRACKING_BUFSIZE, MPI_DOUBLE,
                      masterBuffer, TRACKING_BUFSIZE, MPI_DOUBLE, MPI_COMM_WORLD);
        // each node updates the master beam from the master buffer
        masterBeam->setStepFromBuffer(masterBuffer);
        MPI_Barrier(MPI_COMM_WORLD);
        // now we have all data available on all nodes

        // update the loggers
        if (teufel::rank == 0)
        {
            for (int il=0; il<(int)listLoggers.size(); il++)
                if (listLoggers.at(il)->log_requested(step+1)) listLoggers.at(il)->update();
        }
        
        // handle watch points
        if (teufel::rank == 0)
            for (int iw=0; iw<(int)watches.size(); iw++)
            {
                watch_t w = watches.at(iw);
                if (w.step == step+1)
                {
                    std::cout << "writing watch point " << w.filename << std::endl;
                    int nw = masterBeam->WriteWatchPointHDF5(w.filename);
                    std::cout << nw << " particles written." << std::endl;
                    std::cout << std::endl;
                }
            }
        
        // make a print once every 30s
        current_time = MPI_Wtime();
        if (current_time-print_time > 30)
        {
            if (teufel::rank == 0)
            {
                std::cout << "time elapsed : " << current_time-start_time << " s   ";
                std::cout << "tracking step : " << step << std::endl;
            };
            ss.str(std::string());
            ss << "/proc/" << PID << "/status";
            myfile.open(ss.str());
            if(!myfile.is_open()) {
                cout << "cannot access /proc/{PID}/status - memory usage unknown." << std::endl;
            } else {
                while(getline(myfile, line)) {
                    // std::cout << line << std::endl;
                    size_t pos = line.find(delimiter);
                    std::string name = line.substr(0, pos);
                    std::string rest = line.erase(0,pos+1);
                    if (name=="VmSize") {
                        mem_size = std::stoi(rest);
                    }
                }
                cout << "node " << teufel::rank << " memory usage : " << mem_size/1e3 << " MB" << std::endl;
                myfile.close();
            }
            print_time = current_time;
        };
        
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100000);
    cout << "node " << teufel::rank << " finished tracking " << trackedBeam->getNOP() << " particles"
        << " over " << trackedBeam->getNOTS() << " time steps" << std::endl;
    usleep(100000);
    MPI_Barrier(MPI_COMM_WORLD);
    
    ss.str(std::string());
    ss << "/proc/" << PID << "/status";
    myfile.open(ss.str());
    if(!myfile.is_open()) {
        cout << "cannot access /proc/{PID}/status - memory usage unknown." << std::endl;
    } else {
        while(getline(myfile, line)) {
            // std::cout << line << std::endl;
            size_t pos = line.find(delimiter);
            std::string name = line.substr(0, pos);
            std::string rest = line.erase(0,pos+1);
            if (name=="VmSize") {
                mem_size = std::stoi(rest);
            }
        }
        cout << "node " << teufel::rank << " memory usage : " << mem_size/1e3 << " MB" << std::endl;
        myfile.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100000);
    if (teufel::rank==0) cout << std::endl << std::flush;
    usleep(100000);
    MPI_Barrier(MPI_COMM_WORLD);
    
    // record the finish time
    if (teufel::rank == 0)
    {
        double stop_time = MPI_Wtime();
        std::cout << "time elapsed during tracking : " << stop_time-start_time << " s" << std::endl;
    };

    // write the logged tracking data
    if (teufel::rank == 0)
        for (int i=0; i<(int)listLoggers.size(); i++)
            listLoggers.at(i)->WriteBeamParametersSDDS();

    // ===============================================================
    // Compute the radiation observations :
    // Each node does the computation for it's own set of particles.
    // The fields are combined before writing the files.
    // ===============================================================

    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100000);
    if (teufel::rank==0) cout << std::endl << std::flush;
    usleep(100000);
    MPI_Barrier(MPI_COMM_WORLD);
    
    // compute all observations
    for (int i=0; i<(int)listObservers.size(); i++)
    {
        double start_time = MPI_Wtime();
        if (teufel::rank == 0)
        {
            std::cout << std::endl << "computing observer No. " << i+1 << std::endl;
        };
        Observer *obs = listObservers.at(i);
        MPI_Barrier(MPI_COMM_WORLD);
        usleep(100000);
        std::cout << "node " << teufel::rank << " integrating ... " << std::endl;
        // OpenMP parallelization is handled inside the observer object
        if (obs->getSource()==BeamObservation)
            obs->integrate(trackedBeam);
        if (obs->getSource()==LatticeObservation)
            obs->integrate(lattice);
        double stop_time = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
        usleep(100000);
        if (teufel::rank == 0)
        {
            std::cout << "time elapsed : " << stop_time-start_time << " s" << std::endl;
        };
        
        // collect all the field computed on the individual nodes into the master node
        // TODO: this should be done in chunks, otherwise we temporarily need
        // 3 times the observer memory (original and 2 buffers).
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
        usleep(100000);
        if (teufel::rank == 0)
            std::cout << std::endl << "All buffers allocated." << std::endl << std::flush;
        MPI_Reduce(
            nodeBuffer,                 // send buffer
            reduceBuffer,               // receive buffer
            count,                      // number of values
            MPI_DOUBLE,                 // data type
            MPI_SUM,                    // operation
            0,                          // rank of the root process
            MPI_COMM_WORLD              // communicator
        );
        MPI_Barrier(MPI_COMM_WORLD);
        usleep(100000);
        if (teufel::rank == 0)
            std::cout << "Reduce finished." << std::endl << std::flush;
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
            catch (exception const & e) { cout << e.what() << endl;}
        }
        // if (teufel::rank == 0) delete reduceBuffer;
        delete[] reduceBuffer;
        delete[] nodeBuffer;

        MPI_Barrier(MPI_COMM_WORLD);
        usleep(100000);
        cout << "node " << teufel::rank << " done with observer " << i+1 << std::endl << std::flush;
        usleep(100000);
        MPI_Barrier(MPI_COMM_WORLD);

    }

    // ===============================================================
    // all computation done - clean up
    // ===============================================================

    // delete all observers
    // one has to delete all observers by hand as the list only contains pointer references
    for (int i=0; i<(int)listObservers.size(); i++)
        delete listObservers.at(i);
    listObservers.clear();
    
    // delete all watches
    watches.clear();
    
    // delete all loggers
    // one has to delete all loggers by hand as the list only contains pointer references
    for (int i=0; i<(int)listLoggers.size(); i++)
        delete listLoggers.at(i);
    listLoggers.clear();
    
    // deleting the lattice automatically deletes all lattice elments
    delete lattice;
    
    // deleting the beam automatically deletes all bunches and particles belonging to it
    delete masterBeam;
    delete trackedBeam;

    delete[] masterBuffer;
    delete[] trackingBuffer;

    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100000);
    if (teufel::rank==0)
        std::cout << std::endl << " TEUFEL run finished." << std::endl << std::endl;
    
    MPI_Finalize();
    return 0;
}
