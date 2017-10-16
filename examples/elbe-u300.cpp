/*=========================================================================

  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers

  Module:    undulator test case

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
    \brief Radiation from the ELBE U300 THz source

    @author Ulf Lehnert
    @date 12.4.2017
    @file elbe-u300.cpp
    
    This test case tracks a single particle in an undulator field.
    The particle energy is 24 MeV it corresponds to 4.37e8 electrons
    that is a charge of 70 pC.
    
    That particle propagates through an undulator of 8 periods with 300 mm
    preiod length. The particle starts at z=0, the undulator is centered at z=2.0m.
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
#include "observer.h"
#include "particle.h"
#include "fields.h"
#include "undulator.h"
#include <iostream>
#include <fstream>

int NOP = 1e3;
int NOTS = 4000;    // number of time steps

int main()
{
    double B = 0.100;
    double lambda = 0.300;
    double N = 8;
    PlanarUndulator* Undu = new PlanarUndulator(Vector(0.0, 0.0, 2.0));
    Undu->Setup(B, lambda, N);

    printf("B =  %9.6g T\n", B);
    printf("Undulator Period = %9.6g m\n", lambda);
    printf("N = %9.6g\n", (double)N);
    double gamma = 48;
    double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    double betagamma = sqrt(gamma * gamma - 1.0);
    double K = 0.934 * 100 * lambda * B;
    double lambdar = (lambda / (2 * gamma * gamma)) * (1 + K * K / 2);
    printf("beta =  %12.9g\n", beta);
    printf("gamma =  %12.9g\n", gamma);
    printf("c*p =  %12.9g MeV\n", 1e-6 * mecsquared * betagamma);
    printf("Radiation Wavelength =  %6.3f mm\n", lambdar * 1.0e3);
    printf("Radiation Frequency =  %6.3g THz\n", SpeedOfLight/lambdar*1.0e-12);
    
    // a simple lattice with just the Undulator Field
    Lattice* lattice = new Lattice;
    lattice->addElement(Undu);

    // setup a bunch with 
    // one single particle corresponding to 4.37e8 electrons (70pC)
    // positioned at the origin (by default)
    double ch = 70.0e-12 / ElementaryCharge;
    ChargedParticle *electron = new ChargedParticle(-ch,ch);
    electron->setMomentum(Vector(0.0,0.0,betagamma));
    Bunch *single = new Bunch();
    single->Add(electron);
    
    // an electron bunch of 70pC total charge modeled with NOP particles
    // the initial bunch length is 50fs
    ch = 70.0e-12 / ElementaryCharge / NOP;
    double sigma_t = 50.0e-15;
    double sigma_z = SpeedOfLight * beta * sigma_t;
    printf("sigma_t =  %9.3g ps\n", 1e12*sigma_t);
    printf("sigma_z =  %9.3g mm\n", 1e3*sigma_z);
    Distribution *dist = new Distribution(6, NOP);
    dist->generateGaussian(0.000, 0.001, 0);	// x gaussian with sigma=1mm
    dist->generateGaussian(0.000, 0.0005, 1);	// y gaussian with sigma=0.5mm
    dist->generateGaussian(0.000, sigma_z, 2);	// z gaussian with sigma_z
    dist->generateGaussian(0.000, 0.001, 3);	// px gaussian 1mrad
    dist->generateGaussian(0.000, 0.001, 4);	// py gaussian 1mrad
    dist->generateGaussian(betagamma, 0.001*betagamma, 5);	// pz gaussian 0.1% energy spread
    // set the initial positions and momenta of the particles
    Bunch *bunch = new Bunch(dist, -ch, ch);

    // Tracking should be done for 4.0 m in lab space corresponding to tau [s].
    // Inside the undulator we have an additional pathlength of one radiation
    // wavelength per period. The radiation wavelength already includes the
    // velocity of the particles. Outside the undulator the electron moves with beta*SpeedOfLight
    double tau = (double)N * (lambda + lambdar) / SpeedOfLight + (4.0 - (double)N * lambda) / (beta * SpeedOfLight);
    double deltaT = tau / NOTS;

    // create a particle dump of the inital distribution
    if (0 != bunch->WriteWatchPointSDDS("elbe-u300_start.sdds"))
    {
	printf("SDDS write \033[1;31m failed!\033[0m\n");
    }
    else
    {
	printf("SDDS file written - \033[1;32m OK\033[0m\n");
    }
    
    // track a beam consisting of two bunches
    // one containing only the single (reference) electron
    Beam *beam = new Beam();
    beam->Add(single);
    beam->Add(bunch);
    
    // setup for the tracking procedure
    beam->InitVay(deltaT, lattice);
    // compute the radiation of the single electron, the bunch and the whole beam on axis
    double t0 = 10.0/SpeedOfLight - 1.0e-12;
    PointObserver<Bunch> singleObs = PointObserver<Bunch>(
	single, Vector(0.0, 0.0, 10.0), t0, 0.05e-13, 3000);
    PointObserver<Bunch> bunchObs = PointObserver<Bunch>(
	bunch, Vector(0.0, 0.0, 10.0), t0, 0.05e-13, 3000);
    // PointObserver<Beam> beamObs = PointObserver<Beam>(beam, Vector(0.0, 0.0, 10.0), t0, 0.05e-13, 3000);
    // do the tracking of the beam
    printf("tracking particles ... ");
    fflush(stdout);
    for (int step=0; step<NOTS; step++)
    {
	beam->StepVay(lattice);
	// update the observers every step
	// printf("single ... ");
	// printf(" z(avg.)= %6.3f ", single->avgPosition().z);
	singleObs.integrate();
	// printf("bunch ... ");
	// printf(" z(avg.)= %6.3f ", bunch->avgPosition().z);
	bunchObs.integrate();
	// beamObs->integrate();
	// printf("\n");
    }
    printf("done.\n");
    
    // create a trajectory dump fo the single electron
    /*
    int retval = electron->WriteTrajectorySDDS("elbe-u300_Trajectory.sdds");
    if (0 != retval)
    {
	printf("SDDS write \033[1;31m failed! - error %d\033[0m\n", retval);
    }
    else
    {
        printf("SDDS file written - \033[1;32m OK\033[0m\n");
    }
    */

    // create a particle dump of the final distribution
    if (0 != bunch->WriteWatchPointSDDS("elbe-u300_final.sdds"))
    {
	printf("SDDS write \033[1;31m failed!\033[0m\n");
    }
    else
    {
	printf("SDDS file written - \033[1;32m OK\033[0m\n");
    }
    /*
    if (0 != bunch->WriteWatchPointHDF5(tau, "elbe-u300_final.h5"))
    {
	printf("HDF5 write \033[1;31m failed!\033[0m\n");
    }
    else
    {
	printf("HDF5 file written - \033[1;32m OK\033[0m\n");
    }
    */
    
    // write field time traces
    int retval = singleObs.WriteTimeDomainFieldSDDS("elbe-u300_SingleEl_ObsRadField.sdds");
    if (0 != retval)
    {
	printf("SDDS write \033[1;31m failed! - error %d\033[0m\n", retval);
    }
    else
    {
	printf("SDDS time domain field written - \033[1;32m OK\033[0m\n");
    }
    retval = bunchObs.WriteTimeDomainFieldSDDS("elbe-u300_Bunch_ObsRadField.sdds");
    if (0 != retval)
    {
	printf("SDDS write \033[1;31m failed! - error %d\033[0m\n", retval);
    }
    else
    {
	printf("SDDS time domain field written - \033[1;32m OK\033[0m\n");
    }

    /*
    // write frequency spectrum to file
    retval = Obs.WriteSpectrumSDDS(
	"elbe-u300_RadSpectrum.sdds",
	0.0, 5.0e12, 0.01e12);
    if (0 != retval)
    {
	printf("SDDS spectrum write \033[1;31m failed! - error %d\033[0m\n", retval);
    }
    else
    {
	printf("SDDS spectrum file written - \033[1;32m OK\033[0m\n");
    }
    */
    
    // clean up
    delete lattice;
    // deleting the beam automatically
    // deletes all bunches and particles belonging to it
    delete beam;

    return 0;
}
