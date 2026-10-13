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
    @file teufel_omp.cpp
    
    This is the OMP version of the main executable of TEUFEL.
    It should be executed like
    
    export OMP_NUM_THREADS=8
    teufel_omp input.xml
    
    All computation is distributed over a number of threads on a single node
    in a CPU cluster using OpenMP.
*/

#include <math.h>
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
#include "snapshot.h"
#include "undulator.h"

#include "pugixml.hpp"

int main(int argc, char *argv[])
{

    int NumberOfThreads = omp_get_max_threads();
    
    std::cout << std::endl;
    std::cout << " TEUFEL " << TEUFEL_VERSION_MAJOR << ".";
    std::cout.width(2);
    std::cout.fill('0');
    std::cout << TEUFEL_VERSION_MINOR << ".";
    std::cout.width(2);
    std::cout.fill('0');
    std::cout << TEUFEL_VERSION_PATCH << std::endl;
    cout << std::endl <<" THz-Emission from Undulators and Free-Electron Lasers" << std::endl << std::endl;
    cout << std::endl <<" TEUFEL parallel computing on " << NumberOfThreads << " threads." << std::endl;
    cout << std::endl;
    cout << std::flush;

    char processor_name[256] = {0};
    int mem_size = -1;
    int PID = getpid();
    gethostname(processor_name, sizeof(processor_name)-1);
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

    cout << "host " << processor_name
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
        cout << "memory usage : " << mem_size/1e3 << " MB" << std::endl;
        myfile.close();
    }
    
    cout << std::endl << std::flush;
    
    // ===============================================================
    // parsing the input file
    // ===============================================================
    
    // The first command line argument is interpreted as the input file name.
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
    
    // Further parsing of the input document is done by the parser object.
    InputParser *parse = new InputParser(root);
    
    // We create an empty lattice object.
    // All lattice elements found when parsing the input are added to this.
    Lattice *lattice = new Lattice;
    int NoE = parse->parseLattice(lattice);
    std::cout << std::endl;
    std::cout << "lattice of " << NoE << " elements created." << std::endl;
    cout << std::endl << std::flush;
    
    // We create an empty beam object.
    // Then we call the parser to fill in the necessary information from the input file.
    // The structure of this beam shall then be preserved troughout the tracking procedure.
    // That is necessary to keep the references to the contained bunches for logging
    // or other bunch-related procedures. This way it is possible to log information
    // for inidividual bunches making up the total beam.
    // 
    // The tracked particle information is stored in this object after every
    // tracking step, so the total beam history is available for interactions.
    std::cout << "TEUFEL: creating master beam." << std::endl << std::endl;
    Beam *masterBeam = new Beam();
    std::vector<Logger*> listLoggers;
    ProbeInfo probe{false, 0, "", 0, 0};
    int NoB = parse->parseBeam(masterBeam, &listLoggers, &probe);
    std::cout << std::endl;
    std::cout << "beam of " << NoB << " bunches created." << std::endl;
    std::cout << "total number of particles : " << masterBeam->getNOP() << std::endl;
    std::cout << "total charge : "
        << masterBeam->getTotalCharge()*1.0e9 << "nC" << std::endl;
    std::cout << "defined " << (int)listLoggers.size() << " loggers." << std::endl;
    cout << std::endl << std::flush;

    // get all tracking information from the input file
    std::vector<watch_t> watches;
    std::list<InteractionField*> interactions;
    parse->parseTracking(masterBeam, &watches, &interactions);
    std::cout << std::endl;
    std::cout << "defined " << (int)watches.size() << " watch points." << std::endl;
    // all interaction fields are added to the lattice to be included in the particel tracking
    for (InteractionField* f : interactions) lattice->addElement(f);
    std::cout << "defined " << (int)interactions.size() << " interactions." << std::endl;
    std::cout << "lattice containing " << (int)lattice->count() << " components." << std::endl;
    std::cout << std::endl << std::flush;
    
    // parse all observer definitions
    std::vector<BaseObserver*> listObservers;
    parse->parseObservers(&listObservers);
    std::cout << std::endl;
    std::cout << "defined " << (int)listObservers.size() << " observers." << std::endl;
    std::cout << std::endl << std::flush;

    // We register the number of time steps to track the beam.
    // This has been set by the parser, but will be cleared when
    // initalizing the trajectories for actual tracking.
    int NOTS = masterBeam->getNOTS();
    double time_step = masterBeam->getTimeStep();
    std::cout << "tracking for " << NOTS << " time steps dt=" << time_step << " s" << std::endl;
    std::cout << "master beam trajectory storage : "
        << (double)masterBeam->getNOP() * NOTS * 10.0 * sizeof(double) / 1e6
        << " MB  for " << masterBeam->getNOP() << " particles" << std::endl;
    std::cout << std::endl << std::flush;

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
        cout << "memory usage : " << mem_size/1e3 << " MB" << std::endl;
        myfile.close();
    }
    cout << std::endl << std::flush;

    // ===============================================================
    // reset the master beam to the initial step
    // ===============================================================

    // the initial step is stored in a buffer and re-loaded into the beam
    double *masterBuffer = new double[masterBeam->getStepBufferSize()];
    masterBeam->bufferStep(masterBuffer);

    // update the beam from the buffer
    masterBeam->clearTrajectories();
    // we will use the master beam to store the gathered tracking data
    // so we pre-allocate the necessary trajectory storage for every particle
    masterBeam->preAllocate(NOTS+1);
    masterBeam->setStepFromBuffer(masterBuffer);
    // we have not performed any time step yet
    // all particles just have the step=0 coordinates
    masterBeam->setNOTS(0);
    masterBeam->setTimeStep(time_step);

    std::cout << "master beam setup done." << std::endl;
    std::cout << std::endl << std::flush;

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
        cout << "memory usage : " << mem_size/1e3 << " MB" << std::endl;
        myfile.close();
    }
    cout << std::endl << std::flush;

    // ===============================================================
    // copy all particles into the bunch that will actually be tracked
    // ===============================================================

    // this is the bunch that will actually be tracked
    Bunch *trackedBunch = new Bunch();
    // traverse the beam and copy particles
    int nob = masterBeam->getNOB();
    for (int i=0; i<nob; i++)
    {
        Bunch *b = masterBeam->getBunch(i);
        int nop = b->getNOP();
        for (int j=0; j<nop; j++)
        {
            ChargedParticle *master_particle = b->getParticle(j);
            // create a new particle as an exact copy of the master beam particle
            ChargedParticle *p = new ChargedParticle(master_particle);
            p->preAllocate(NOTS+1);
            trackedBunch->Add(p);
        }
    }
    
    // ===============================================================
    // track the beam
    // ===============================================================

    // this is the beam we will actually track
    Beam *trackedBeam = new Beam();
    trackedBeam->Add(trackedBunch);
    // copy the tracking information from the master beam
    trackedBeam->setTrackingMethod(masterBeam->getTrackingMethod());
    trackedBeam->setTimeStep(masterBeam->getTimeStep());
    // prepare the tracking of the beam
    trackedBeam->setupTracking(lattice);

    // setup the ProbeLogger to observe trackedBeam
    if (probe.requested)
    {
        std::cout << "writing " << probe.number << " particles ";
        std::cout << "every " << probe.step << " steps ";
        std::cout << "to " << probe.filename << std::endl;
        std::cout << std::endl;
        ProbeLogger<Beam>* pl = new ProbeLogger<Beam>(trackedBeam, probe.filename, probe.step, probe.number);
        listLoggers.push_back(pl);
    };
    
    std::cout << "TEUFEL: tracking a beam of " <<
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
        cout << "memory usage : " << mem_size/1e3 << " MB" << std::endl;
        myfile.close();
    }
    
    cout << std::endl << std::flush;

    // handle watch point of initial particle distribution if requested
    for (int iw=0; iw<(int)watches.size(); iw++)
    {
        watch_t w = watches.at(iw);
        if (w.step == 0)
        {
            std::cout << "writing watch point " << w.filename << std::endl;
            if (w.type == hdf5)
            {
                int nw = masterBeam->WriteWatchPointHDF5(w.filename);
                std::cout << nw << " particles written." << std::endl;
                std::cout << std::endl;
            };
            if (w.type == sdds)
            {
                int res = masterBeam->WriteWatchPointSDDS(w.filename);
                if (res==0)
                {
                    std::cout << masterBeam->getNOP() << " particles written." << std::endl;
                    std::cout << std::endl;
                } else {
                    std::cout << "error " << res << " while writing." << std::endl;
                    std::cout << std::endl;
                };
            };
        }
    }

    // update loggers with initial particle distribution
    for (int il=0; il<(int)listLoggers.size(); il++)
        listLoggers.at(il)->update();
                    
    // initalize the interaction fields
    // field source is the master beam - it is queried before every tracking step
    // the fields act on the tracked beam
    std::cout << "TEUFEL: initializing interactions." << std::endl;
    for (InteractionField* f : interactions) f->init(masterBeam);
    
    // record the start time
    double start_time = omp_get_wtime();
    double print_time = start_time;
    double current_time = start_time;

    // this is the global tracking clock
    // if the beam internal clock deviates this will result in exceptions from the interaction fields
    double tracking_time = 0.0;
    double tracking_time_step = masterBeam->getTimeStep();
    std::cout << "start tracking : beam is at step=" << masterBeam->getNOTS();
    std::cout << " t=" << tracking_time << " s" << std::endl;
    std::cout << std::endl;

    // do the tracking of the beam
    // Shared memory parallelization is used for field integration only.
    // It may be employed in the future during tracking for
    // computing the interaction field between the particles.
    // With a-priori fields tracking is not parallelized.
    for (int step=0; step<NOTS; step++)
    {
    
        // compute the interactions
        // source is the master beam containing all particles
        //! @todo we should not compute the interaction from the master beam
        //! but from the tracked beam
        for (InteractionField* f : interactions) f->update(tracking_time);
        
        // do a step
        trackedBeam->doStep(lattice);
        tracking_time += tracking_time_step;
        
        // for (int i=0; i<trackedBunch->getNOP(); i++)
        // {
        //     ChargedParticle *p = trackedBunch->getParticle(i);
        //     p->StepVay(lattice);
        // }
            
        // the master beam is updated from the tracked beam
        trackedBeam->bufferStep(masterBuffer);
        masterBeam->setStepFromBuffer(masterBuffer);
        // stepping sanity check
        if (masterBeam->getNOTS()!=step+1)
        {
            std::cout << "  beam is at step=" << masterBeam->getNOTS() << "  main() is at step=" << step << std::endl;
            throw(IOexception("main() - beam is at wrong step"));
        }

        // update the loggers
        for (int il=0; il<(int)listLoggers.size(); il++)
            if (listLoggers.at(il)->log_requested(step+1)) listLoggers.at(il)->update();
        
        // handle watch points
        for (int iw=0; iw<(int)watches.size(); iw++)
        {
            watch_t w = watches.at(iw);
            if (w.step == step+1)
            {
                std::cout << "writing watch point " << w.filename << std::endl;
                if (w.type == hdf5)
                {
                    int nw = masterBeam->WriteWatchPointHDF5(w.filename);
                    std::cout << nw << " particles written." << std::endl;
                    std::cout << std::endl;
                };
                if (w.type == sdds)
                {
                    int res = masterBeam->WriteWatchPointSDDS(w.filename);
                    if (res==0)
                    {
                        std::cout << masterBeam->getNOP() << " particles written." << std::endl;
                        std::cout << std::endl;
                    } else {
                        std::cout << "error " << res << " while writing." << std::endl;
                        std::cout << std::endl;
                    };
                };
            }
        }
        
        // make a print once every 30s
        current_time = omp_get_wtime();
        if (current_time-print_time > 30)
        {
            std::cout << "time elapsed : " << current_time-start_time << " s   ";
            std::cout << "tracking step : " << step << std::endl;
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
                cout << "memory usage : " << mem_size/1e3 << " MB" << std::endl;
                myfile.close();
            }
            print_time = current_time;
        };
        
    }
    
    cout << "finished tracking " << trackedBeam->getNOP() << " particles"
        << " over " << trackedBeam->getNOTS() << " time steps" << std::endl;
    
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
        cout << "memory usage : " << mem_size/1e3 << " MB" << std::endl;
        myfile.close();
    }
    cout << std::endl << std::flush;
    
    // record the finish time
    double stop_time = omp_get_wtime();
    std::cout << "time elapsed during tracking : " << stop_time-start_time << " s" << std::endl << std::endl;

    // write the logged data
    for (int i=0; i<(int)listLoggers.size(); i++)
        listLoggers.at(i)->WriteData();

    for (InteractionField* f : interactions)
        f->write_output();

    // ===============================================================
    // Compute the radiation observations
    // ===============================================================

    cout << std::endl << std::flush;
    
    // we remove all interaction fields from the lattice
    // we don't want these to be included with the lattice observations
    int i_lat = 0;
    while (i_lat<lattice->count())
    {
        GeneralField* elem = lattice->at(i_lat);
        if (!elem->is_external())
        {
            std::cout << "deleting interaction field from lattice at index " << i_lat  << std::endl;
            lattice->remove(i_lat);
        } else { 
            i_lat++;
        }
    }
    
    //! @todo one should not delete the interaction fields but move them
    //! to an own list. Then an observer could look at either the beam,
    //! the lattice or the interaction fields
    
    // compute all observations
    for (int i=0; i<(int)listObservers.size(); i++)
    {
        double start_time = omp_get_wtime();
        std::cout << std::endl << "computing observer No. " << i+1 << std::endl;
        BaseObserver *observer = listObservers.at(i);
        std::cout << "integrating ... " << std::endl;
        // OpenMP parallelization is handled inside the observer object
        // try and handle the different observer types
        if (auto obs = dynamic_cast<Observer<Beam>*>(observer))
        {
            obs->setSource(trackedBeam);
            obs->integrate();
        }
        if (auto obs = dynamic_cast<Observer<Lattice>*>(observer))
        {
            obs->setSource(lattice);
            obs->integrate();
        }
        double stop_time = omp_get_wtime();
        std::cout << "time elapsed : " << stop_time-start_time << " s" << std::endl;
        
        // a SnapshotObserver need special attention - it needs to
        // buffer all particle positions from the master beam
        if (SnapshotObserver<Beam> *obs = dynamic_cast<SnapshotObserver<Beam>*>(observer))
        {
            obs->setSource(masterBeam);
            std::cout << "SnaphotObserver collecting particle coordinates from master beam." << std::endl << std::flush;
            int num = obs->storeParticlePositions();
            std::cout << "SnaphotObserver stored coordinates of " << num << " particles." << std::endl << std::flush;
        }
        // produce the output
        try
        { 
            observer->generateOutput();
            std::cout << "Observer output file written." << std::endl;
        }
        catch (exception const & e) { cout << e.what() << endl;}

        cout << "done with observer " << i+1 << std::endl << std::flush;

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

    std::cout << std::endl << " TEUFEL run finished." << std::endl << std::endl;
    
    return 0;
}
