/*=========================================================================

  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers

  Module:    undulator example (single particle)

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

/*!
    \brief Radiation emitted by a point charge in a waveguided undulator

    @author Ulf Lehnert
    @date 27.8.2019
    @file waveguide_wiggler.cpp
    
    This test case tracks a number of particles in an undulator field.
    
    That particles propagate through an undulator of 40 periods with 110 mm
    preiod length. The particle starts at z=0, the undulator is centered at z=3.0m.
    At z=6m the produced radiation is observed. A parallel plate waveguide is simulated
    by mirror particles of the tracked particle. With a sufficient number
    of mirror trajectories a waveguide up to the observation screen is obtained.
    
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
#include "global.h"
#include "logger.h"
#include "observer.h"
#include "particle.h"
#include "screen.h"
#include "fields.h"
#include "undulator.h"

int NOP = 100;      // number of particles
int NOTS = 6000;	// number of time steps for tracking
int NOM = 50;       // number of mirror reflections
int NX = 51;        // screen size in horizontal direction
double dX = 0.005;  // screen pixel size [m] in horizontal direction
int NY = 101;       // screen size in vertical direction
double dY = 0.0001; // screen pixel size [m] in vertical direction
int NOSTS = 6000;	// number of time steps for the screen

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
        std::cout << " TEUFEL " << std::endl;
        cout << std::endl <<" THz-Emission From Undulators and Free-Electron Lasers" << std::endl << std::endl;
        cout << std::endl <<" TEUFEL parallel computing on " << NumberOfCores << " cores." << std::endl << std::endl;
    }

    // ===============================================================
    // setup and tracking is done by all processes in parallel
    // only the node rank 0 will print output
    // ===============================================================
    
    double gamma = 14.0/0.511;
    double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    double betagamma = sqrt(gamma * gamma - 1.0);
    if (teufel::rank==0) printf("beta =  %12.9g\n", beta);
    if (teufel::rank==0) printf("gamma =  %12.9g\n", gamma);
    if (teufel::rank==0) printf("c*p =  %12.9g MeV\n", 1e-6 * mecsquared * betagamma);

    double lambdar = 0.01/85.0;
    if (teufel::rank==0) printf("Radiation Wavelength =  %6.3f µm\n", lambdar * 1.0e6);
    if (teufel::rank==0) printf("Radiation Wavenumber =  %6.3f 1/cm\n", 0.01/lambdar);
    if (teufel::rank==0) printf("Radiation Frequency =  %6.3g THz\n", SpeedOfLight/lambdar*1.0e-12);

    double lambda = 0.110;
    double N = 40;
    if (teufel::rank==0) printf("Undulator Period = %9.6g m\n", lambda);
    if (teufel::rank==0) printf("N = %9.6g\n", (double)N);
    double Ksq = lambdar/lambda * 2*gamma*gamma;
    double K = sqrt(2.0*(Ksq-1.0));
    double B = K / 93.4 / lambda;
    if (teufel::rank==0) printf("B =  %9.6g T\n", B);
    PlanarUndulator* Undu = new PlanarUndulator(Vector(0.0, 0.0, 3.0));
    Undu->Setup(B, lambda, N);

    // a simple lattice with just the Undulator Field
    Lattice* lattice = new Lattice;
    lattice->addElement(Undu);

    // setup a bunch with 1nC charge
    double ch = 1.0e-9 / ElementaryCharge;
    
    /*
    // one single particle of 1nC
    // positioned at the origin (by default)
    ChargedParticle *electron = new ChargedParticle(-ch,ch);
    electron->initTrajectory(0.0,
        Vector(0.0,0.0,0.0),
        Vector(0.0,0.0,betagamma),
        Vector(0.0,0.0,0.0));
    Bunch *primary = new Bunch();
    primary->Add(electron);
    */

    // setup a bunch with many particles but zero bunchlength
    // match transverse beam properties to undulator for 50 mm mrad normalized emittance
    Distribution *dist = new Distribution(6, NOP);
    // set the initial positions and momenta of the particles
    double lambda_b = gamma*sqrt(2.0)/K * lambda;
    double sig_y = sqrt(50e-6/betagamma * lambda_b/(2.0*M_PI));
    double sig_yp = 50e-6/betagamma / sig_y;
    if (teufel::rank==0)
    {
        printf("Betatron Period = %9.6g m\n", lambda_b);
        printf("matched beam size sigma_y = %f mm\n",sig_y*1000.0);
        printf("matched beam momentum sigma_yp = %f mrad\n",sig_yp*1000.0);
    }
    // transverse emittance 1.83µm (geometric) 50µm (normalized)
    dist->generateGaussian(0, 0.0);	    // x gaussian
    dist->generateGaussian(1, sig_y);	// y gaussian
    dist->generateGaussian(2, 0.0);	    // z gaussian
    dist->generateGaussian(3, 0.0);	    // px gaussian
    dist->generateGaussian(4, sig_yp*betagamma);	// py gaussian
    dist->generateGaussian(5, 0.0);	    // pz gaussian
    // particles are "back transported" by 0.8m (undulator entrance to start)
    // by adding a -0.8 m/rad correlated position change (4 acts on 1)
    dist->addCorrelation(4, 1, -0.8/betagamma);
    Bunch *primary = new Bunch(dist, 0.0, Vector(0.0,0.0,0.0), Vector(0.0,0.0,betagamma), -ch, ch);
    
    // Tracking should be done for 5.5 m in lab space corresponding to tau [s].
    // Inside the undulator we have an additional pathlength of one radiation
    // wavelength per period. The radiation wavelength already includes the
    // velocity of the particles. Outside the undulator the electron moves with beta*SpeedOfLight
    double tau = (double)N * (lambda + lambdar) / SpeedOfLight + (5.5 - (double)N * lambda) / (beta * SpeedOfLight);
    double deltaT = tau / NOTS;

    // setup for the tracking procedure
    primary->InitVay(deltaT, lattice);

    // create a particle dump of the initial distribution
    if (teufel::rank==0)
    {
        if (0 != primary->WriteWatchPointSDDS("Waveguide_initial.sdds"))
        {
	        printf("SDDS write \033[1;31m failed!\033[0m\n");
        }
        else
        {
        	printf("SDDS file written - \033[1;32m OK\033[0m\n");
        }
    };
    // log the Parameters of the bunch
    TrackingLogger<Bunch> *bunchLog = new TrackingLogger<Bunch>(primary, "Waveguide_BeamParam.sdds");

    // do the tracking of the beam
    if (teufel::rank==0) printf("tracking particles ...\n");
    if (teufel::rank==0) fflush(stdout);
    // record the start time
    timespec start_time, stop_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);
    for (int step=0; step<NOTS; step++)
    {
    	primary->StepVay(lattice);
    	// log every 10th step
	    if (step % 10 == 0) bunchLog->update();
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
	
    // record the finish time
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop_time);
    if (teufel::rank==0) printf("finished tracking particles.\n");
    double elapsed = stop_time.tv_sec-start_time.tv_sec + 1e-9*(stop_time.tv_nsec-start_time.tv_nsec);
    if (teufel::rank==0) printf("time elapsed during tracking : %9.3f s\n",elapsed);

    // create a particle dump of the final distribution
    if (teufel::rank==0)
    {
        if (0 != primary->WriteWatchPointSDDS("Waveguide_final.sdds"))
        {
	        printf("SDDS write \033[1;31m failed!\033[0m\n");
        }
        else
        {
        	printf("SDDS file written - \033[1;32m OK\033[0m\n");
        }
    };

    // create a tracking parameter dump fo the bunch
    if (teufel::rank==0)
    {
        int retval = bunchLog->WriteBeamParametersSDDS("Waveguide_BeamParam.sdds");
        if (0 != retval)
        {
	        printf("SDDS write \033[1;31m failed! - error %d\033[0m\n", retval);
        }
        else
        {
            printf("SDDS file written - \033[1;32m OK\033[0m\n");
        }
    };

    // ===============================================================
    // The beam of the node rank 0 is distributed to all nodes
    // so the nodes can do independent computations on an
    // identical set of particle trajectories
    // ===============================================================

    // buffers for MPI communication
    // we rely on all particles having equal trajectory length
    int nop = primary->getNOP();
    ChargedParticle *p = primary->getParticle(0);
    int trajsize = p->getNP();
    int bufsize = p->TrajBufSize();
    double *particlebuffer = new double[bufsize];
    if (teufel::rank == 0) printf("broadcasting %d particles\n",nop);

    // the particles are copied from the primary bunch into this one
    Bunch *trackedBunch = new Bunch();
    // all nodes synchronously loop over the beam
    for (int j=0; j<nop; j++)
    {
        p = primary->getParticle(j);
        // in principle only root would have to buffer its particle trajectory
        p->serializeTraj(particlebuffer);
        // broadcast the buffer from root to all nodes
        MPI_Bcast(particlebuffer, bufsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        ChargedParticle *copy = new ChargedParticle(particlebuffer, trajsize);
        trackedBunch->Add(copy);
    }
    if (teufel::rank == 0) std::cout << "broadcast finished." << std::endl;
    delete particlebuffer;
    
    // ===============================================================
    // Compute radiation emitted by the tracked beam distributed over all nodes.
    // After the computation the fields are gatered onto the root node.
    // Only the root node then writes the output file.
    // ===============================================================
    
	// compute the radiation observed on a finite screen 3.0m downstream the undulator center
	// this observer only sees the direct particle
	double t0 = (6.0 - 10*lambdar) / SpeedOfLight;
	double dt = lambdar / 12.0 / SpeedOfLight;
    ScreenObserver *directObs = new ScreenObserver(
        "Waveguide_MPI_direct.h5",
    	Vector(0.0, 0.0, 6.0),      // position
    	Vector(dX, 0.0, 0.0),       // dx
    	Vector(0.0, dY, 0.0),       // dy
    	NX,                         // unsigned int nx,
    	NY,                         // unsigned int ny,
    	t0,
    	dt,
    	NOSTS);                     // number of time steps
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);
    directObs->integrate_mp(trackedBunch, NumberOfCores, teufel::rank);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop_time);
    if (teufel::rank==0)
    {
        printf("finished parallel radiation calculation on screen.\n");
        elapsed = stop_time.tv_sec-start_time.tv_sec + 1e-9*(stop_time.tv_nsec-start_time.tv_nsec);
        printf("time elapsed : %9.3f s\n",elapsed);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    // collect all the field computed on the individual nodes into the root node
    unsigned int count = directObs->getBufferSize();
    std::cout << "Node " << teufel::rank << " allocating buffers for "<< count << " doubles" << std::endl;
    // fill the buffer and get its address
    double* nodeBuffer = directObs->getBuffer();
    if (nodeBuffer==0)
    {
        std::cout << "MPI_Reduce for Obsserver : node " << teufel::rank << " was unable to get the node buffer." << std::endl;
        throw(IOexception("memory allocation error"));
    };
    // we have to define a buffer for the sum on all nodes
    double* reduceBuffer = new double[count];
    if (reduceBuffer==0)
    {
        std::cout << "MPI_Reduce for Obsserver : node " << teufel::rank << " was unable to get the reduce buffer." << std::endl;
        throw(IOexception("memory allocation error"));
    };
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "node " << teufel::rank << " : nodeBuffer " << nodeBuffer << " reduceBuffer " << reduceBuffer << std::endl;
    if (teufel::rank == 0)
        std::cout << std::endl << "All buffers allocated." << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
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
    // the root node writes the output file    
    if (teufel::rank==0)
    {
        std::cout << "Reduce finished." << std::endl;
        // the root node copies the data from the reduce buffer into the Observer
        directObs->fromBuffer(reduceBuffer,count);
        // write screen time traces
        try
        { 
        	directObs->WriteTimeDomainFieldHDF5();
        	printf("Screen observer time domain field written - \033[1;32m OK\033[0m\n");
        }
        catch (exception& e) { cout << e.what() << endl;}
    };
    
    delete reduceBuffer;
    delete nodeBuffer;
    delete directObs;
    MPI_Barrier(MPI_COMM_WORLD);

    // ===============================================================
    // Create a bunch containing all mirror particles needed to simulate
    // waveguided emission in addition to the tracked particles.
    // ===============================================================

    MPI_Barrier(MPI_COMM_WORLD);
    if (teufel::rank==0) printf("computing mirror particles.\n");
    nop = trackedBunch->getNOP();
    Bunch *mirrorBunch = new Bunch();
    for (int j=0; j<nop; j++)
    {
        // copy the primary particle
        ChargedParticle *p = trackedBunch->getParticle(j);
        ChargedParticle *e0 = new ChargedParticle(p);
        mirrorBunch->Add(e0);
    	// the waveguide of 10mm height extends up to infinity
    	// mirror planes are at y= +/- 5mm, 10mm, 15mm, ... with reversed polarity
    	// even numbers of mirrors constitute a shift without polarity change
    	for (int m=1; m<=NOM; m++)
    	{
            ChargedParticle *e1 = new ChargedParticle(p);
            ChargedParticle *e2 = new ChargedParticle(p);
            if (m%2 == 0)
            {
                e1->Shift(Vector(0.0,0.010*m,0.0), 1.0);
                e2->Shift(Vector(0.0,-0.010*m,0.0), 1.0);
            }
            else
            {
                e1->Mirror(Vector(0.0,0.005*m,0.0), Vector(0.0,1.0,0.0), -1.0);
                e2->Mirror(Vector(0.0,-0.005*m,0.0), Vector(0.0,-1.0,0.0), -1.0);
            };
            mirrorBunch->Add(e1);
            mirrorBunch->Add(e2);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // ===============================================================
    // Compute radiation emitted by the mirror beam distributed over all nodes.
    // After the computation the fields are gatered onto the root node.
    // Only the root node then writes the output file.
    // ===============================================================

    if (teufel::rank==0) printf("computing radiation on screen including mirror particles.\n");
	// compute the radiation observed on a finite screen 3.0m downstream the undulator center
	// this observer only sees the direct particle
    ScreenObserver *mirrorObs = new ScreenObserver(
        "Waveguide_MPI_mirror.h5",
    	Vector(0.0, 0.0, 6.0),      // position
    	Vector(dX, 0.0, 0.0),       // dx
    	Vector(0.0, dY, 0.0),       // dy
    	NX,                         // unsigned int nx,
    	NY,                         // unsigned int ny,
    	t0,
    	dt,
    	NOSTS);                     // number of time steps
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);
    mirrorObs->integrate_mp(mirrorBunch, NumberOfCores, teufel::rank);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop_time);
    if (teufel::rank==0)
    {
        printf("finished parallel radiation calculation on screen.\n");
        elapsed = stop_time.tv_sec-start_time.tv_sec + 1e-9*(stop_time.tv_nsec-start_time.tv_nsec);
        printf("time elapsed : %9.3f s\n",elapsed);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    
    // collect all the field computed on the individual nodes into the root node
    count = mirrorObs->getBufferSize();
    std::cout << "Node " << teufel::rank << " allocating buffers for "<< count << " doubles" << std::endl;
    // fill the buffer and get its address
    nodeBuffer = mirrorObs->getBuffer();
    if (nodeBuffer==0)
    {
        std::cout << "MPI_Reduce for Obsserver : node " << teufel::rank << " was unable to get the node buffer." << std::endl;
        throw(IOexception("memory allocation error"));
    };
    // we have to define a buffer for the sum on all nodes
    reduceBuffer = new double[count];
    if (reduceBuffer==0)
    {
        std::cout << "MPI_Reduce for Obsserver : node " << teufel::rank << " was unable to get the recude buffer." << std::endl;
        throw(IOexception("memory allocation error"));
    };
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "node " << teufel::rank << " : nodeBuffer " << nodeBuffer << " reduceBuffer " << reduceBuffer << std::endl;
    if (teufel::rank == 0)
        std::cout << std::endl << "All buffers allocated." << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
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
    // the root node writes the output file    
    if (teufel::rank==0)
    {
        std::cout << "Reduce finished." << std::endl;
        // the root node copies the data from the reduce buffer into the Observer
        mirrorObs->fromBuffer(reduceBuffer,count);
        // write screen time traces
        try
        { 
        	mirrorObs->WriteTimeDomainFieldHDF5();
        	printf("Screen observer time domain field written - \033[1;32m OK\033[0m\n");
        }
        catch (exception& e) { cout << e.what() << endl;}
    };
    
    delete reduceBuffer;
    delete nodeBuffer;
    delete mirrorObs;
    MPI_Barrier(MPI_COMM_WORLD);

    // clean up
    delete lattice;
    // deleting the beam automatically
    // deletes all bunches and particles belonging to it
    delete primary;
    delete trackedBunch;
    delete mirrorBunch;
   
    if (teufel::rank==0)
        std::cout << std::endl << " TEUFEL run finished." << std::endl << std::endl;
    
    MPI_Finalize();
    return 0;
}
