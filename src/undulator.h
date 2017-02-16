/*=========================================================================
 * 
 *  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers
 * 
 *  Copyright (c) 2017 U. Lehnert
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * =========================================================================*/

#ifndef UNDULATOR_H
#define UNDULATOR_H

#include "vector.h"
#include "externalfield.h"

using namespace std;

// a class for a planar undulator

// the beam direction is along the z axis starting at z=0
// the main magnetic field component is in y direction
// so the undulator deflects in x direction
// the field profile is flat in x direction (no focusssing)

class Undulator : public ExternalField
{

  public:

    // constructor
    Undulator(double B,                         // peak field [T]
              double lambda,                    // undulator period [m]
              int    N                          // number of undulator periods
             );

    double  GetBPeak();
    double  GetLambdaU();
    int     GetNPeriods();
    double  GetKpeak();
    double  GetKrms();

  private:

    Vector ElementLocalEField(double t, Vector X);
    Vector ElementLocalBField(double t, Vector X);

    double  BPeak;                              // peak field [T]
    double  LambdaU;                            // undulator period [m]
    int     NPeriods;                           // number of undulator periods
    double  Krms;                               // undulator parameter
    double  ky, kz;                             // undulator periodicity [1/m]
};

#endif