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
    
    The program generates a trajectory dump elbe-u300_trajectory.sdds which
    can be used to plot the elctron trajectory.
    
    At c*t=2.0m a snapshot of the particle distribution and the local fields is generated.

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>

#include "bunch.h"
#include "beam.h"
#include "global.h"
#include "logger.h"
#include "observer.h"
#include "particle.h"
#include "fields.h"
#include "undulator.h"

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <time.h>

int NOP = 1e3;		// number of particles
int NOTS = 4000;	// number of time steps

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);
    
    int NumberOfCores = 1;
    int my_rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &NumberOfCores);
    bool parallel = (NumberOfCores != 1);
    if (parallel) MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if (parallel and my_rank==0)
    {
        printf("\n  TEUFEL parallel computing on %d cores.\n\n", NumberOfCores);
    }

    double B = 0.25242;
    double lambda = 0.300;
    double N = 8;
    PlanarUndulator* Undu = new PlanarUndulator(Vector(0.0, 0.0, 2.0));
    Undu->Setup(B, lambda, N);

    if (my_rank==0)
    {
        printf("Undulator Period = %9.6g m\n", lambda);
        printf("N = %9.6g\n", (double)N);
        printf("B =  %9.6g T\n", B);
        printf("K(rms) =  %9.3g\n", Undu->GetKrms());
    }
    double gamma = 195.695;
    double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    double betagamma = sqrt(gamma * gamma - 1.0);
    double K = Undu->GetKpeak();
    double lambdar = (lambda / (2 * gamma * gamma)) * (1 + K * K / 2);
    if (my_rank==0)
    {
        printf("beta =  %12.9g\n", beta);
        printf("gamma =  %12.9g\n", gamma);
        printf("c*p =  %12.9g MeV\n", 1e-6 * mecsquared * betagamma);
        printf("Radiation Wavelength =  %6.3f mm\n", lambdar * 1.0e3);
        printf("Radiation Frequency =  %6.3g THz\n", SpeedOfLight/lambdar*1.0e-12);
    }
        
    // a simple lattice with just the Undulator Field
    Lattice* lattice = new Lattice;
    lattice->addElement(Undu);

    // an electron bunch of 1nC total charge modeled with NOP particles
    // the initial bunch length is 100fs
    double ch = 451.0e-12 / ElementaryCharge / NOP;
    double sigma_t = 0.0450e-12;
    double sigma_z = SpeedOfLight * beta * sigma_t;
    Distribution *dist = new Distribution(6, NOP);
    // allocate a communication buffer for the distribution
    int bufsize = 6*NOP;
    double *buffer = new double[bufsize];
    // only the master node computes the particle distribution
    if (my_rank==0)
    {
        printf("sigma_t =  %9.3g ps\n", 1e12*sigma_t);
        printf("sigma_z =  %9.3g mm\n", 1e3*sigma_z);
        printf("\n");
        // set the initial positions and momenta of the particles
        // transverse emittance 0.051µm (geometric) 10µm (normalized)
        dist->generateGaussian(0.000, 0.000309, 0);	// x gaussian with sigma=1.0mm
        dist->generateGaussian(0.000, 0.000165*betagamma, 3);	// px gaussian px/pz=0.131mrad
        // particles are "back transported" by 2m (undulator center to start)
        // by adding a -2 m/rad correlated position change
        dist->addCorrelation(3, 0, -2.0/betagamma);
        dist->generateGaussian(0.000, 0.000309, 1);	// y gaussian with sigma=0.7mm
        dist->generateGaussian(0.000, 0.000165*betagamma, 4);	// py gaussian py/pz=0.131mrad
        // particles are "back transported" by 0.8m (undulator entrance to start)
        // by adding a -0.8 m/rad correlated position change
        dist->addCorrelation(4, 1, -0.8/betagamma);
        dist->generateGaussian(0.000, sigma_z, 2);	// z gaussian with sigma_z
        dist->generateGaussian(betagamma, 0.001*betagamma, 5);	// pz gaussian 0.1% energy spread
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

    // node 1 generates a bunch of all particles for output only
    /*
    if (my_rank==1)
    {
        Bunch *all = new Bunch(dist, -ch, ch);
        if (0 != all->WriteWatchPointSDDS("MPI_node1_all_particles_received.sdds"))
        {
          	printf("SDDS write \033[1;31m failed!\033[0m\n");
        }
        else
        {
    	    printf("SDDS file written - \033[1;32m OK\033[0m\n");
        }
        delete all;
    }
    */
    
    // The complete set of particles is divided into an number of bunches.
    // Every bunch is tracked on one computation node all running in parallel.
    int bunchsize = NOP / NumberOfCores;
    int remainder = NOP % NumberOfCores;   
    int index;  // where do my particles start in the list
    if (my_rank<remainder)
    {   // first few nodes get one particle more
        bunchsize++;
        index = my_rank*(bunchsize+1);
    }
    else 
    {   
        index = remainder + my_rank*bunchsize;
    }
    Distribution *mydist = dist->subDist(index, bunchsize);
    printf("Node #%d computing %d particles starting at index %d\n",my_rank,bunchsize,index);
    Bunch *bunch = new Bunch(mydist, -ch, ch);

    // every node creates a dump o the initial particle distribution
    /*
	char partFileName[80];
	sprintf(partFileName,"MPI_node%d_particles_start.sdds",my_rank);
    if (0 != bunch->WriteWatchPointSDDS(partFileName))
    {
      	printf("SDDS write \033[1;31m failed!\033[0m\n");
    }
    else
    {
	    printf("SDDS file written - \033[1;32m OK\033[0m\n");
    }
    */
    
    // we are done with the distributions
    delete buffer;
    delete dist;
    delete mydist;
    
    // Tracking should be done for 4.0 m in lab space corresponding to tau [s].
    // Inside the undulator we have an additional pathlength of one radiation
    // wavelength per period. The radiation wavelength already includes the
    // velocity of the particles. Outside the undulator the electron moves with beta*SpeedOfLight
    double tau = (double)N * (lambda + lambdar) / SpeedOfLight + (4.0 - (double)N * lambda) / (beta * SpeedOfLight);
    double deltaT = tau / NOTS;

    // track a beam containing just the one bunch per node
    Beam *beam = new Beam();
    beam->Add(bunch);
    
    // setup for the tracking procedure
    beam->InitVay(deltaT, lattice);
    // Record the radiation of the beam at 10m distance from the undulator center.
    // Every node just records the radiation emitted by its own particles.
    // We have to sum it up later.
    double z0 = 2.0 + 10.0;
    double t0 = z0/SpeedOfLight - 1.0e-12;
    ScreenObserver<Bunch> screenObs = ScreenObserver<Bunch>(
    	bunch,
    	Vector(0.0, 0.0, z0),		// position
    	Vector(0.007, 0.0, 0.0),		// dx
    	Vector(0.0, 0.007, 0.0),		// dy
    	81,				// unsigned int nx,
    	81,				// unsigned int ny,
    	t0,
    	2.0e-14,			// double dt,
    	2000);				// NOTS

    // log the Parameters of the bunch
    /*
    TrackingLogger<Bunch> *bunchLog = new TrackingLogger<Bunch>(bunch);
    */
    
    // do the tracking of the beam
    printf("Node #%d tracking started.\n",my_rank);
    // record the start time
    double start_time = MPI_Wtime();
    double print_time = start_time;
    for (int step=0; step<NOTS; step++)
    {
    	beam->StepVay(lattice);
	    // update the observers every step
	    screenObs.integrate();
	    // log every 10th step
	    // if (step % 10 == 0) bunchLog->update();
	    // make a print once every 60s
	    double current_time = MPI_Wtime();
	    if (current_time-print_time > 120)
	    {
	        print_time = current_time;
	        printf("Node #%d tracking step %d / %d\n",my_rank,step,NOTS);
	    };
    }
    // record the finish time
    double stop_time = MPI_Wtime();
    printf("Node #%d finished after %6.2f s.\n",my_rank,stop_time-start_time);
	
	// create a file for the field per node
	/*
	char obsFileName[80];
	sprintf(obsFileName,"MPI_node%d_Screen_ObsRadField.h5",my_rank);
	screenObs.WriteTimeDomainFieldHDF5(obsFileName);
	*/
	
	/*
    // create a tracking parameter dump fo the bunch
    int retval = bunchLog->WriteBeamParametersSDDS("dali-u300_BeamParam.sdds");
    if (0 != retval)
    {
    	printf("SDDS write \033[1;31m failed! - error %d\033[0m\n", retval);
    }
    else
    {
        printf("SDDS file written - \033[1;32m OK\033[0m\n");
    }

    // create a particle dump of the final distribution
    if (0 != bunch->WriteWatchPointSDDS("dali-u300_final.sdds"))
    {
    	printf("SDDS write \033[1;31m failed!\033[0m\n");
    }
    else
    {
	    printf("SDDS file written - \033[1;32m OK\033[0m\n");
    }
    */
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank == 0) printf("\nAll nodes finished tracking.\n\n");
    // collect all the field computed on the individual nodes
    // into the master node
    unsigned int count = screenObs.getCount();
    printf("Node #%d allocating buffers for %d doubles\n",my_rank,count);
    // fill the buffer and get its address
    double* nodeBuffer = screenObs.getBuffer();
    if (nodeBuffer==0)
        printf("MPI_Reduce for screenObs : node #%d was unable to get the node buffer.\n",my_rank);
    // we have to define a buffer for the sum on all nodes
    double* reduceBuffer = new double[count];
    if (reduceBuffer==0)
        printf("MPI_Reduce for screenObs : node #%d was unable to get the reduce buffer.\n",my_rank);
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank == 0) printf("\nAll buffers allocated.\n");
    MPI_Reduce(
        nodeBuffer,                 // send buffer
        reduceBuffer,               // receive buffer
        count,                      // number of values
        MPI_DOUBLE,                 // data type
        MPI_SUM,                    // operation
        0,                          // rank of the root process
        MPI_COMM_WORLD              // communicator
    );
    if (my_rank == 0) printf("\nReduce finished.\n\n");
    // the root node copies the data from the reduce buffer into the sreenObs
    // and writes it to the output file
    if (my_rank == 0)
    {
        screenObs.fromBuffer(reduceBuffer,count);
        try
        { 
        	screenObs.WriteTimeDomainFieldHDF5("MPI_dali-u300_Screen_ObsRadField.h5");
        	printf("Screen observer time domain field written - \033[1;32m OK\033[0m\n");
        }
        catch (exception& e) { cout << e.what() << endl;}
    }
    // if (my_rank == 0) delete reduceBuffer;
    delete reduceBuffer;
    delete nodeBuffer;
        
    // clean up
    delete lattice;
    // deleting the beam automatically deletes all bunches and particles belonging to it
    delete beam;
    // delete observers and loggers
    // delete bunchLog;

    MPI_Barrier(MPI_COMM_WORLD);
    if (parallel and my_rank==0)
    {
        printf("\n  TEUFEL run finished.\n\n");
    }
    
    MPI_Finalize();
    return 0;
}
