/*=========================================================================

  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers

  Module:    Main Program

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "particle.h"
#include "externalfield.h"
#include "homogeneousmagnet.h"

#include "SDDS.h"

// <bunch> parameters
int NOP = 1;                    // number of (real) particles
int NOTS = 3000;                // number of time steps

int main ()
{
  FILE *dump;
  int d;                // integer buffer used for file writing

  printf("\nTEUFEL Version 0.01 U.Lehnert 12/2016\n");
  printf("homogeneous magnet testcase\n\n");

  Vector B=Vector(0.033166247903554,0.05,0.08);
  HomogeneousMagnet *mag = new HomogeneousMagnet(B);
  printf("B =  %9.6g T\n",(mag->getB0()).norm());
  printf("Bx = %9.6g T  ",(mag->getB0()).x);
  printf("By = %9.6g T  ",(mag->getB0()).y);
  printf("Bz = %9.6g T\n",(mag->getB0()).z);
  
  double cp = 10.0e6;   // 10 MeV
  printf("c*p =  %9.6g MeV\n",1e-6*cp);
  double betagamma= cp / mecsquared;
  printf("gamma =  %9.6g\n",sqrt(1.0+betagamma*betagamma));
  double Rgyro = betagamma * mecsquared / SpeedOfLight / mag->getB0().norm();
  printf("R =  %9.6g m\n",Rgyro);

  Lattice *lattice = new Lattice;
  lattice->addElement(mag);

  ChargedParticle *bunch = new ChargedParticle[NOP];

  for (int i=0; i<NOP; i++)
    {
      // initial position of the particle
      // amplitude of fig-8 oscillation is a0/k
      Vector X0 = Vector(0.0, 0.0, Rgyro);
      // initial momentum of the particle
      Vector P0 = Vector(betagamma, 0.0, 0.0);
      // bunch[i].TrackLF(NOTS, 5.0e-12, X0, P0, lattice);
      bunch[i].TrackVay(NOTS, 5.0e-12, X0, P0, lattice);
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

  // write trajectories to SDDS file
  SDDS_DATASET SDDS_dataset;

  /* SDDS_InitializeOutput
      arguments:
       1) *SDDS_DATASET
       2) data mode = SDDS_ASCII or SDDS_BINARY
       3) lines per row of ascii data
       4) description string
       5) contents string
       6) output file name
      return:
       1 on success
       0 on failure 
  */
  if (SDDS_InitializeOutput(&SDDS_dataset, SDDS_BINARY, 1,
                            NULL, NULL, "traj.sdds" ) != 1) {
    return(1);
  }
  fprintf(stdout, "output initialized\n");

  /* SDDS_DefineSimpleParameter
      arguments:
       1) *SDDS_DATASET
       2) parameter name
       3) parameter units
       4) SDDS datatype
      return:
       1 on success
       0 on failure 
  */
  if (SDDS_DefineSimpleParameter(&SDDS_dataset, "B", "T", SDDS_DOUBLE)!=1 ||
      SDDS_DefineSimpleParameter(&SDDS_dataset, "cp",  "eV", SDDS_DOUBLE)!=1 ||
      SDDS_DefineSimpleParameter(&SDDS_dataset, "gamma", NULL, SDDS_DOUBLE)!=1 ||
      SDDS_DefineSimpleParameter(&SDDS_dataset, "R", "m", SDDS_DOUBLE)!=1) {
    return(1);
  }
  fprintf(stdout, "parameters defined\n");

  /* SDDS_WriteLayout
      arguments:
       1) *SDDS_DATASET
      return:
       1 on success
       0 on failure 
  */
  if (SDDS_WriteLayout(&SDDS_dataset)!=1) {
    return(1);
  }
  fprintf(stdout, "layout written\n");

  /* SDDS_StartPage
      arguments:
       1) *SDDS_DATASET
       2) expected number of rows
      return:
       1 on success
       0 on failure 
  */
  if (SDDS_StartPage(&SDDS_dataset, 5)!=1) {
    return(1);
  }
  fprintf(stdout, "started page\n");

  /* SDDS_SetParameters
      arguments:
       1) *SDDS_DATASET
       2) mode = SDDS_SET_BY_NAME | SDDS_PASS_BY_VALUE
        3) char *name1, value1, char *name2, value2, ..., NULL
       2) mode = SDDS_SET_BY_NAME | SDDS_PASS_BY_REFERENCE
        3) char *name1, void *data1, char *name2, void *data2, ..., NULL
       2) mode = SDDS_SET_BY_INDEX | SDDS_PASS_BY_VALUE
        3) long index1, value1, long index2, value2, ..., -1
       2) mode = SDDS_SET_BY_INDEX | SDDS_PASS_BY_REFERENCE
        3) long index1, void *data1, long index2, void *data2, ..., -1
      return:
       1 on success
       0 on failure 
  */
  if (SDDS_SetParameters(&SDDS_dataset, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, 
			 "B", mag->getB0().norm(), 
                         "cp",  cp, 
                         "gamma", sqrt(1.0+betagamma*betagamma), 
                         "R",  Rgyro,
                         NULL)!=1) {
    return(1);
  }
  fprintf(stdout, "parameters set\n");

  /* SDDS_WritePage
      arguments:
       1) *SDDS_DATASET
      return:
       1 on success
       0 on failure 
  */
  if (SDDS_WritePage(&SDDS_dataset)!=1) {
    return(1);
  }
  fprintf(stdout,"page written out\n");

  /* SDDS_Terminate
      arguments:
       1) *SDDS_DATASET
      return:
       1 on success
       0 on failure 
  */
  if (SDDS_Terminate(&SDDS_dataset)!=1) {
    return(1);
  }

}
