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
    teufel::my_rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &NumberOfCores);
    MPI_Comm_rank(MPI_COMM_WORLD, &teufel::my_rank);
    
    if (teufel::my_rank==0)
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
    
    // The first command line argument is interpreted as the input file name
    // This is an XML document which is opened and parsed here.
    pugi::xml_document doc;
    if (argc < 2) {
        if (teufel::my_rank==0) std::cout << "Usage is teufel <infile>" << std::endl << std::endl;
        exit(1);
    } else {
        if (teufel::my_rank==0) std::cout << "reading XML input from " << argv[1] << std::endl;
        pugi::xml_parse_result result = doc.load_file(argv[1]);
        if (result)
        {
            if (teufel::my_rank==0) std::cout << "input parsed without errors" << std::endl;
        }
        else
        {
            if (teufel::my_rank==0) std::cout << "ERROR reading file " << argv[1] << std::endl;
            if (teufel::my_rank==0) std::cout << "Error description: " << result.description() << std::endl;
            exit(1);
        }
    };
    pugi::xml_node root = doc.child("teufel");
    if (!root)
        throw(IOexception("TEUFEL::InputParser - root node <teufel> not found."));
    string description = root.attribute("description").value();
    string author = root.attribute("author").value();
    if (teufel::my_rank==0) std::cout << "case : " << description << std::endl;
    if (teufel::my_rank==0) std::cout << "by : " << author << std::endl << std::endl;
    
    // Further parsing of the input document is done by the simulation object
    InputParser *parse = new InputParser(root);
    
    // We create an empty lattice object.
    // All lattice elements found when parsing the input are added to this.
    Lattice *lattice = new Lattice;
    int NoE = parse->parseLattice(lattice);
    if (teufel::my_rank==0) std::cout << std::endl;
    if (teufel::my_rank==0) std::cout << "lattice of " << NoE << " elements created." << std::endl;
    
    // We create an empty beam object.
    // Then we call the parser to fill in the necessary information
    // from the input file.
    Beam *beam = new Beam();
    int NoB = parse->parseBeam(beam);
    if (teufel::my_rank==0) std::cout << std::endl;
    if (teufel::my_rank==0) std::cout << "beam of " << NoB << " bunches created." << std::endl;
    if (teufel::my_rank==0) std::cout << "total number of particles : " << beam->getNOP() << std::endl;
    if (teufel::my_rank==0) std::cout << "total charge : " << beam->getTotalCharge()*ElementaryCharge*1.0e9 << "nC" << std::endl;
    if (teufel::my_rank==0) std::cout << std::endl;

    // get all tracking information from the input file
    std::vector<watch_t> watches;
    int NoW = parse->parseTracking(beam, &watches);
    if (teufel::my_rank==0) std::cout << "defined " << NoW << " watch points." << std::endl;
    if (NoW != (int)watches.size())
        if (teufel::my_rank==0) std::cout << "WARNING : number of watches (" << NoW << ") differs from length of the list (" << watches.size() << ")" << std::endl;
    if (teufel::my_rank==0) std::cout << std::endl;
    
    // parse all observer definitions
    std::vector<Observer*> listObservers;
    int NoO = parse->parseObservers(&listObservers, beam);
    if (teufel::my_rank==0) std::cout << std::endl;
    if (teufel::my_rank==0) std::cout << "defined " << NoO << " observers." << std::endl;
    if (NoO != (int)listObservers.size())
        if (teufel::my_rank==0) std::cout << "WARNING : number of observers (" << NoO << ") differs from length of the list (" << listObservers.size() << ")" << std::endl;
    if (teufel::my_rank==0) std::cout << std::endl;
    
    // we are done with the input document
    delete parse;

    // ===============================================================
    // done parsing the input file
    // ===============================================================

/*
    double B = 0.10315;
    double lambda = 0.300;
    double N = 8;
    PlanarUndulator* Undu = new PlanarUndulator(Vector(0.0, 0.0, 2.0));
    Undu->Setup(B, lambda, N);

    if (teufel::my_rank==0)
    {
        printf("Undulator Period = %9.6g m\n", lambda);
        printf("N = %9.6g\n", (double)N);
        printf("B =  %9.6g T\n", B);
        printf("K(rms) =  %9.3g\n", Undu->GetKrms());
    }
    double gamma = 50.881;
    double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    double betagamma = sqrt(gamma * gamma - 1.0);
    double K = Undu->GetKpeak();
    double lambdar = (lambda / (2 * gamma * gamma)) * (1 + K * K / 2);
    if (teufel::my_rank==0)
    {
        printf("beta =  %12.9g\n", beta);
        printf("gamma =  %12.9g\n", gamma);
        printf("c*p =  %12.9g MeV\n", 1e-6 * mecsquared * betagamma);
        printf("Radiation Wavelength =  %6.3f mm\n", lambdar * 1.0e3);
        printf("Radiation Frequency =  %6.3g THz\n", SpeedOfLight/lambdar*1.0e-12);
    }
        
    // a simple lattice with just the Undulator Field
    // Lattice* lattice = new Lattice;
    // lattice->addElement(Undu);

    // an electron bunch of 77pC total charge modeled with NOP particles
    // the initial bunch length is 100fs
    double ch = 77.0e-12 / ElementaryCharge / NOP;
    double sigma_t = 0.102e-12;
    double sigma_z = SpeedOfLight * beta * sigma_t;
    Distribution *dist = new Distribution(6, NOP);
    // allocate a communication buffer for the distribution
    int bufsize = 6*NOP;
    double *buffer = new double[bufsize];
    // only the master node computes the particle distribution
    if (teufel::my_rank==0)
    {
        printf("sigma_t =  %9.3g ps\n", 1e12*sigma_t);
        printf("sigma_z =  %9.3g mm\n", 1e3*sigma_z);
        printf("\n");
        // set the initial positions and momenta of the particles
        // transverse emittance 0.051µm (geometric) 10µm (normalized)
        dist->generateGaussian(0.000, 0.000483, 0);     // x gaussian with sigma=1.0mm
        dist->generateGaussian(0.000, 0.000407*betagamma, 3);   // px gaussian px/pz=0.131mrad
        // particles are "back transported" by 2m (undulator center to start)
        // by adding a -2 m/rad correlated position change
        dist->addCorrelation(3, 0, -2.0/betagamma);
        dist->generateGaussian(0.000, 0.000483, 1);     // y gaussian with sigma=0.7mm
        dist->generateGaussian(0.000, 0.000407*betagamma, 4);   // py gaussian py/pz=0.131mrad
        // particles are "back transported" by 0.8m (undulator entrance to start)
        // by adding a -0.8 m/rad correlated position change
        dist->addCorrelation(4, 1, -0.8/betagamma);
        dist->generateGaussian(0.000, sigma_z, 2);      // z gaussian with sigma_z
        dist->generateGaussian(betagamma, 0.001*betagamma, 5);  // pz gaussian 0.1% energy spread
        // the master node knows the distribution and fills its buffer
        dist->bufferData(buffer,bufsize);
        // generate a bunch of all particles for output only
        Bunch *all = new Bunch(dist, -ch, ch);
        if (0 != all->WriteWatchPointSDDS("MPI_node0_all_particles_generated.sdds"))
        {
                printf("SDDS write \033[1;31m failed!\033[0m\n");
        }
        else
        {
            printf("SDDS file written - \033[1;32m OK\033[0m\n");
        }
        delete all;
    }

    // Distribute the particle coordinates to all nodes.
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(buffer, bufsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    // Now all nodes should have the same buffer.
    dist->fromBuffer(buffer,bufsize);

    // The complete set of particles is divided into an number of bunches.
    // Every bunch is tracked on one computation node all running in parallel.
    int bunchsize = NOP / NumberOfCores;
    int remainder = NOP % NumberOfCores;   
    int index;  // where do my particles start in the list
    if (teufel::my_rank<remainder)
    {   // first few nodes get one particle more
        bunchsize++;
        index = teufel::my_rank*(bunchsize+1);
    }
    else 
    {   
        index = remainder + teufel::my_rank*bunchsize;
    }
    Distribution *mydist = dist->subDist(index, bunchsize);
    printf("Node #%d computing %d particles starting at index %d\n",teufel::my_rank,bunchsize,index);
    Bunch *bunch = new Bunch(mydist, -ch, ch);

    // we are done with the distributions
    delete buffer;
    delete dist;
    delete mydist;
    
    // Tracking should be done for 3.4 m in lab space corresponding to tau [s].
    // Inside the undulator we have an additional pathlength of one radiation
    // wavelength per period. The radiation wavelength already includes the
    // velocity of the particles. Outside the undulator the electron moves with beta*SpeedOfLight
    double tau = (double)N * (lambda + lambdar) / SpeedOfLight + (3.4 - (double)N * lambda) / (beta * SpeedOfLight);
    double deltaT = tau / NOTS;

    // track a beam containing just the one bunch per node
    // Beam *beam = new Beam();
    beam->Add(bunch);
    
    // setup for the tracking procedure
    beam->setTimeStep(deltaT);
    beam->InitVay(lattice);

    // log the Parameters of the bunch
    
    // do the tracking of the beam
    printf("Node #%d tracking started.\n",teufel::my_rank);
    // record the start time
    double start_time = MPI_Wtime();
    double print_time = start_time;
    for (int step=0; step<NOTS; step++)
    {
        beam->StepVay(lattice);
            // log every 10th step
            // if (step % 10 == 0) bunchLog->update();
            // make a print once every 60s
            double current_time = MPI_Wtime();
            if (current_time-print_time > 120)
            {
                print_time = current_time;
                printf("Node #%d tracking step %d / %d\n",teufel::my_rank,step,NOTS);
            };
    }
    // record the finish time
    double stop_time = MPI_Wtime();
    printf("Node #%d finished after %6.2f s.\n",teufel::my_rank,stop_time-start_time);
        
    // create a file for the field per node
        
    MPI_Barrier(MPI_COMM_WORLD);
    if (teufel::my_rank == 0) printf("\nAll nodes finished tracking.\n\n");
    
    // Record the radiation of the beam at 1.625m distance from the undulator center.
    // Every node just records the radiation emitted by its own particles.
    // We have to sum it up later.
    double z0 = 2.0 + 1.625;
    double t0 = z0/SpeedOfLight - 1.0e-12;
    ScreenObserver<Bunch> screenObs = ScreenObserver<Bunch>(
        bunch,
        "MPI_elbe-u300_Screen_ObsRadField.h5",
        Vector(0.0, 0.0, z0),           // position
        Vector(0.002, 0.0, 0.0),                // dx
        Vector(0.0, 0.002, 0.0),                // dy
        41,                             // unsigned int nx,
        41,                             // unsigned int ny,
        t0,
        5.0e-13,                        // double dt,
        50);                            // NOTS

    // compute field seen from the bunch
    // in parallel on every node for its own particles
    printf("Node #%d integrating...\n",teufel::my_rank);
    start_time = MPI_Wtime();
    screenObs.integrate();
    // record the finish time
    stop_time = MPI_Wtime();
    printf("Node #%d finished after %6.2f s.\n",teufel::my_rank,stop_time-start_time);

    MPI_Barrier(MPI_COMM_WORLD);
    if (teufel::my_rank == 0) printf("\nAll nodes finished the field computation.\n\n");
    
    // collect all the field computed on the individual nodes
    // into the master node
    unsigned int count = screenObs.getCount();
    printf("Node #%d allocating buffers for %d doubles\n",teufel::my_rank,count);
    // fill the buffer and get its address
    double* nodeBuffer = screenObs.getBuffer();
    if (nodeBuffer==0)
        printf("MPI_Reduce for screenObs : node #%d was unable to get the node buffer.\n",teufel::my_rank);
    // we have to define a buffer for the sum on all nodes
    double* reduceBuffer = new double[count];
    if (reduceBuffer==0)
        printf("MPI_Reduce for screenObs : node #%d was unable to get the reduce buffer.\n",teufel::my_rank);
    MPI_Barrier(MPI_COMM_WORLD);
    if (teufel::my_rank == 0) printf("\nAll buffers allocated.\n");
    MPI_Reduce(
        nodeBuffer,                 // send buffer
        reduceBuffer,               // receive buffer
        count,                      // number of values
        MPI_DOUBLE,                 // data type
        MPI_SUM,                    // operation
        0,                          // rank of the root process
        MPI_COMM_WORLD              // communicator
    );
    if (teufel::my_rank == 0) printf("\nReduce finished.\n\n");
    // the root node copies the data from the reduce buffer into the sreenObs
    // and writes it to the output file
    if (teufel::my_rank == 0)
    {
        screenObs.fromBuffer(reduceBuffer,count);
        try
        { 
                screenObs.WriteTimeDomainFieldHDF5();
                printf("Screen observer time domain field written - \033[1;32m OK\033[0m\n");
        }
        catch (exception& e) { cout << e.what() << endl;}
    }
    // if (teufel::my_rank == 0) delete reduceBuffer;
    delete reduceBuffer;
    delete nodeBuffer;
*/

    // clean up
    delete lattice;
    // deleting the beam automatically deletes all bunches and particles belonging to it
    delete beam;
    // delete observers and loggers
    // delete bunchLog;

    MPI_Barrier(MPI_COMM_WORLD);
    if (teufel::my_rank==0)
    {
        printf("\n  TEUFEL run finished.\n\n");
    }
    
    MPI_Finalize();
    return 0;
}
