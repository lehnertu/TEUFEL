/*=========================================================================

  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers

  Module:    Diffraction screen test case

  Copyright (c) 2020 U. Lehnert

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
    \brief Fields of a free electron recorded on an annular diffraction screen

    @author Ulf Lehnert
    @date 3.4.2020
    @file teufel.DiffractionScreen.cpp
    
    This case tests the handling of meshed screens. The file describing such
    a screen is provided as DiffractionScreen.h5 in the tests/ folder.
    It is generated from a Jupyter notebook (Python) CreateDiffractionScreen.ipynb
    provided in the scripts/ folder.
    
    This test case tracks a single electron from the origin until it has crossed
    the plane of a diffraction screen at (0,0,1)m. The screen is angled 45 deg
    to the incoming beam as to reflect radiation into the x direction.
    The geometry is read from the file and after recording the fields
    all data is written back to the same file (which increases in size).
    
    The test is considered passed if there are no errors recorded
    during the file handling.
    
    @return The number of errors encountered in the above list of checks is reported.
    
 */

#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "mesh.h"
#include "particle.h"
#include "fields.h"

// number of time steps
#define NOTS 20000

// the test is supposed to be run by the script in the teufel root directory
// the files to be used
#define INFILE "tests/DiffractionScreen.h5"
#define OUTFILE "tests/DiffractionScreenWithFields.h5"

int main ()
{

    printf("\nTEUFEL - diffraction screen testcase\n");

    // charge [C] of the particle
    double charge = -100e-12;

    // energy of the particle
    double gamma = 20.0;

    // an empty lattice
    Lattice *lattice = new Lattice;
    
    // one bunch containing a single charge representing many electrons
    ChargedParticle *electrons = new ChargedParticle(charge/ElementaryCharge,charge/ElementaryCharge);
    Bunch *bunch = new Bunch();
    bunch->Add(electrons);
    
    // initial position at the origin
    Vector X0 = Vector(0.0, 0.0, 0.0);
    // initial momentum of the particle
    double beta = sqrt(1.0-1.0/(gamma*gamma));
    double betagamma = sqrt(gamma * gamma - 1.0);
    Vector P0 = Vector(0.0, 0.0, betagamma);
    Vector A0 = Vector(0.0, 0.0, 0.0);
    electrons->initTrajectory(0.0, X0, P0, A0);
    
    // count the errors
    int errors = 0;

    // track for 2 m
    double duration = 2.0/(beta*SpeedOfLight);
    double dt = duration / NOTS;

    // track the particle
    bunch->InitVay(dt, lattice);
    for (int i=0; i<NOTS; i++) bunch->StepVay(lattice);
    double t_final = electrons->getTime();
    Vector XP = electrons->getPosition();
    printf("x(%9.6g s) =  (%9.6g,%9.6g,%9.6g) m\n",t_final,XP.x,XP.y,XP.z);

    // check final particle position
    if (fabs(XP.z-2.0) > 1e-6) {
		errors++;
		printf("error in final position %9.6g m - \033[1;31m test failed!\033[0m\n", XP.z);
    };
    
    // record fields on the screen
    try
    {
        // read the screen setup from file
        MeshedScreen* screen = new MeshedScreen(INFILE);
        screen->init();
        screen->zero();
        screen->writeReport(&cout);
        // accumulate the fields
        // TODO: set fields to zero before computation,
        // otherwise we will add to the possibly already existing fields
        printf("computing fields ...\n");
        screen->integrate(bunch);
        // write data
        screen->writeReport(&cout);
        screen->writeFile(OUTFILE);
        // clean up
        delete screen;
    }
    catch (exception& e)
    {
        errors++;
        cout << e.what() << endl;
    }
    
    // clean up
    delete lattice;
    delete electrons;

    return errors;
}
