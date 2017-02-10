/*=========================================================================

  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers

  Module:    Homogeneous magnetic field test case

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
    \brief Homogeneous magnetic field test case

    This test case tracks a single electron in a homogeneous magnetic field.
    This electron moves on an periodic circular trajectory.
    The cyclotron frequency and trajectory radius are compared to known values.

    @author Ulf Lehnert
    @date 10.2.2017
    @file teufel.magnet.cpp
 
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "particle.h"
#include "externalfield.h"

// <bunch> parameters
int NOP = 1;                    // number of (real) particles
int NOTS = 3000;                // number of time steps

class HomogeneousMagnet : public ExternalField
{

  public:

    // constructor
    HomogeneousMagnet(double B)                 // peak field [T]
    { BPeak = B; };

  private:

    Vector ElementLocalEField(double t, Vector X)
    { return Vector(0.0 ,0.0 ,0.0); };

    Vector ElementLocalBField(double t, Vector X)
    { return Vector(0.0 ,BPeak ,0.0); };

  public:

    double  BPeak;                              // peak field [T]
};


int main ()
{
  FILE *dump;
  int d;                // integer buffer used for file writing

  printf("\nTEUFEL Version 0.01 U.Lehnert 10/2014\n");
  printf("homogeneous magnet testcase\n\n");

  // define oscillator
  double B = 0.4;
  HomogeneousMagnet *mag = new HomogeneousMagnet(B);
  printf("B =  %9.6g T\n",mag->BPeak);
  double cp = 10.0e6;   // 10 MeV
  printf("c*p =  %9.6g MeV\n",1e-6*cp);
  double betagamma= cp / mecsquared;
  printf("gamma =  %9.6g\n",sqrt(1.0+betagamma*betagamma));
  double Rgyro = betagamma * mecsquared / SpeedOfLight / B;
  printf("R =  %9.6g m\n",Rgyro);

  ChargedParticle *bunch = new ChargedParticle[NOP];

  for (int i=0; i<NOP; i++)
    {
      // initial position of the particle
      // amplitude of fig-8 oscillation is a0/k
      Vector X0 = Vector(0.0, 0.0, Rgyro);
      // initial momentum of the particle
      Vector P0 = Vector(betagamma, 0.0, 0.0);
      bunch[i].TrackLF(NOTS, 1.0e-12, X0, P0, mag);
    }

  // dump the particle trajectories into a binary file
  dump = fopen("particles.bin","wb");
  // first write the number of particles
  d = NOP;
  fwrite(&d, sizeof(int), 1, dump);
  // write the number of time steps
  d = NOTS;
  fwrite(&d, sizeof(int), 1, dump);
  // write the traces
  // it is transposed while writing because z (the 1-step index)
  // ist used as the outer loop index
  for (int ip=0; ip<NOP; ip++)
    for (int it=0; it<NOTS; it++)
      {
        Vector XP = bunch[ip].TrajPoint(it);
        // fwrite expects a pointer to the value to be written
        fwrite(&(XP.x), sizeof(double), 1, dump);
        fwrite(&(XP.y), sizeof(double), 1, dump);
        fwrite(&(XP.z), sizeof(double), 1, dump);
         // write the momentum in addition
        XP = bunch[ip].TrajMomentum(it);
        fwrite(&(XP.x), sizeof(double), 1, dump);
        fwrite(&(XP.y), sizeof(double), 1, dump);
        fwrite(&(XP.z), sizeof(double), 1, dump);
     };
  fclose(dump);


}
