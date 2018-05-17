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

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "config.h"
#include "simulation.h"

#include "pugixml.hpp"

int main(int argc, char* argv[])
{
    std::cout << std::endl;
    std::cout << " TEUFEL " << TEUFEL_VERSION_MAJOR << ".";
    std::cout.width(2);
    std::cout.fill('0');
    std::cout << TEUFEL_VERSION_MINOR << ".";
    std::cout.width(2);
    std::cout.fill('0');
    std::cout << TEUFEL_VERSION_PATCH << std::endl;
    // printf("\n TEUFEL %d.%02d.%02d\n", TEUFEL_VERSION_MAJOR,TEUFEL_VERSION_MINOR,TEUFEL_VERSION_PATCH);
    cout << std::endl <<" THz-Emission From Undulators and Free-Electron Lasers" << std::endl << std::endl;
    
    // The first command line argument is interpreted as the input file name
    // This is an XML document which is opened and parsed here.
    pugi::xml_document doc;
    if (argc < 2) {
        std::cout << " Usage is teufel <infile>\n\n";
        exit(1);
    } else {
        std::cout << " reading XML input from " << argv[1] << std::endl;
        pugi::xml_parse_result result = doc.load_file(argv[1]);
        if (result)
        {
            std::cout << " input parsed without errors" << std::endl;
        }
        else
        {
            std::cout << " ERROR reading file " << argv[1] << std::endl;
            std::cout << " Error description: " << result.description() << std::endl;
            exit(1);
        }
    };
    pugi::xml_node root = doc.child("teufel");
    if (!root) throw std::invalid_argument("root node <teufel> not found");
    string description = root.attribute("description").value();
    string author = root.attribute("author").value();
    std::cout << " case : " << description << std::endl;
    std::cout << " by : " << author << std::endl << std::endl;
    
    // Further parsing of the input document is done by the simulation object
    Simulation *sim = new Simulation(root);
    int NoE = sim->parseLattice();
    
/*    
    double B = 0.384;
    double lambda = 0.300;
    double N = 8;
    PlanarUndulator* Undu = new PlanarUndulator(Vector(0.0, 0.0, 2.0));
    Undu->Setup(B, lambda, N);

    printf("Undulator Period = %9.6g m\n", lambda);
    printf("N = %9.6g\n", (double)N);
    printf("B =  %9.6g T\n", B);
    printf("K(rms) =  %9.3g\n", Undu->GetKrms());
    double gamma = 195.695;
    double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    double betagamma = sqrt(gamma * gamma - 1.0);
    double K = Undu->GetKpeak();
    double lambdar = (lambda / (2 * gamma * gamma)) * (1 + K * K / 2);
    printf("beta =  %12.9g\n", beta);
    printf("gamma =  %12.9g\n", gamma);
    printf("c*p =  %12.9g MeV\n", 1e-6 * mecsquared * betagamma);
    printf("Radiation Wavelength =  %6.3f mm\n", lambdar * 1.0e3);
    printf("Radiation Frequency =  %6.3g THz\n", SpeedOfLight/lambdar*1.0e-12);
    
    // a simple lattice with just the Undulator Field
    Lattice* lattice = new Lattice;
    lattice->addElement(Undu);

    // an electron bunch of 1nC total charge modeled with NOP particles
    // the initial bunch length is 100fs
    double ch = 1000.0e-12 / ElementaryCharge / NOP;
    double sigma_t = 0.100e-12;
    double sigma_z = SpeedOfLight * beta * sigma_t;
    printf("sigma_t =  %9.3g ps\n", 1e12*sigma_t);
    printf("sigma_z =  %9.3g mm\n", 1e3*sigma_z);
    Distribution *dist = new Distribution(6, NOP);
    // set the initial positions and momenta of the particles
    // transverse emittance 0.051µm (geometric) 10µm (normalized)
    dist->generateGaussian(0.000, 0.000251, 0); // x gaussian with sigma=1.0mm
    dist->generateGaussian(0.000, 0.000204*betagamma, 3);       // px gaussian px/pz=0.131mrad
    // particles are "back transported" by 2m (undulator center to start)
    dist->addCorrelation(3, 0, -2.0/betagamma);
    dist->generateGaussian(0.000, 0.000251, 1); // y gaussian with sigma=0.7mm
    dist->generateGaussian(0.000, 0.000204*betagamma, 4);       // py gaussian py/pz=0.131mrad
    // particles are "back transported" by 0.8m (undulator entrance to start)
    dist->addCorrelation(4, 1, -0.8/betagamma);
    dist->generateGaussian(0.000, sigma_z, 2);  // z gaussian with sigma_z
    dist->generateGaussian(betagamma, 0.001*betagamma, 5);      // pz gaussian 0.1% energy spread
    Bunch *bunch = new Bunch(dist, -ch, ch);
    delete dist;

    // Tracking should be done for 4.0 m in lab space corresponding to tau [s].
    // Inside the undulator we have an additional pathlength of one radiation
    // wavelength per period. The radiation wavelength already includes the
    // velocity of the particles. Outside the undulator the electron moves with beta*SpeedOfLight
    double tau = (double)N * (lambda + lambdar) / SpeedOfLight + (4.0 - (double)N * lambda) / (beta * SpeedOfLight);
    double deltaT = tau / NOTS;

    // create a particle dump of the inital distribution
    if (0 != bunch->WriteWatchPointSDDS("dali-u300_start.sdds"))
    {
      	printf("SDDS write \033[1;31m failed!\033[0m\n");
    }
    else
    {
	    printf("SDDS file written - \033[1;32m OK\033[0m\n");
    }
    
    // track a beam
    Beam *beam = new Beam();
    beam->Add(bunch);
    
    // setup for the tracking procedure
    beam->InitVay(deltaT, lattice);
    // record the radiation of the beam at 10m distance from the undulator center
    double z0 = 2.0 + 10.0;
    double t0 = z0/SpeedOfLight - 1.0e-12;
    PointObserver<Bunch> bunchObs = PointObserver<Bunch>(
	bunch, Vector(0.0, 0.0, z0), t0, 0.05e-13, 4000);
    ScreenObserver<Bunch> screenObs = ScreenObserver<Bunch>(
	bunch,
	Vector(0.0, 0.0, z0),		// position
	Vector(0.01, 0.0, 0.0),		// dx
	Vector(0.0, 0.01, 0.0),		// dy
	81,				// unsigned int nx,
	81,				// unsigned int ny,
	t0,
	0.2e-13,			// double dt,
	2000);				// NOTS

    // log the Parameters of the bunch
    TrackingLogger<Bunch> *bunchLog = new TrackingLogger<Bunch>(bunch);
    
    // do the tracking of the beam
    printf("tracking particles ...\n");
    fflush(stdout);
    // record the start time
    timespec start_time, current_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);
    timespec print_time = start_time;
    for (int step=0; step<NOTS; step++)
    {
    	beam->StepVay(lattice);
	    // update the observers every step
	    // printf("bunch ... ");
	    // printf(" z(avg.)= %6.3f ", bunch->avgPosition().z);
	    bunchObs.integrate();
	    screenObs.integrate();
	    // printf("\n");
	    // log every 10th step
	    if (step % 10 == 0) bunchLog->update();
	    // make a print once every 10s
	    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &current_time);
	    double elapsed = current_time.tv_sec-print_time.tv_sec +
	        1e-9*(current_time.tv_nsec-print_time.tv_nsec);
            if (elapsed>10.0)
	    {
	        print_time = current_time;
	        printf("tracking step %d / %d\n",step,NOTS);
	    };
    }
    // record the finish time
    timespec stop_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop_time);
    printf("finished tracking particles.\n");
    double elapsed = stop_time.tv_sec-start_time.tv_sec +
	1e-9*(stop_time.tv_nsec-start_time.tv_nsec);
    printf("time elapsed during tracking : %9.3f s\n",elapsed);
	
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
    
    // write field time traces
    try
    {
    	bunchObs.WriteTimeDomainFieldSDDS("dali-u300_Bunch_ObsRadField.sdds");
    	printf("SDDS time domain field written - \033[1;32m OK\033[0m\n");
    }
    catch (exception& e) { cout << e.what() << endl;}
    try
    { 
    	screenObs.WriteTimeDomainFieldHDF5("dali-u300_Screen_ObsRadField.h5");
    	printf("Screen observer time domain field written - \033[1;32m OK\033[0m\n");
    }
    catch (exception& e) { cout << e.what() << endl;}

    // clean up
    delete lattice;
    // deleting the beam automatically deletes all bunches and particles belonging to it
    delete beam;
    // delete observers and loggers
    delete bunchLog;
*/

    delete sim;
    
    return 0;
}
