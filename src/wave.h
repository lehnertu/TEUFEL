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

#ifndef WAVE_H
#define WAVE_H

#include "vector.h"
#include "externalfield.h"

using namespace std;

// a class for an electromagnetic plane wave in vacuum
// the wave moves along the z axis
// the polarization direction is along the x axis

class PlaneWave : public ExternalField
{

  public:

    // constructor
    PlaneWave( double lambda,           // wavelength [m]
               double intensity );      // power density [W/cm²]

    double Omega();                     // report the frequency [Hz]
    double EPeak();                     // report the maximum electric field strength [V/m]
    double BPeak();                     // report the maximum magnetic field strength [T]
    double A0();                        // report the dimensionless intensity

  private:

    double Lambda;                              // wavelength [m]
    double omega;                               // frequency [Hz]
    double Intensity;                           // power density (averaged pointing vector) [W/m²]
    double Epeak;                               // peak electric field [V/m]
    double Bpeak;                               // peak magnetic field [T]
    double a0;                                  // dimensionless intensity
    Vector k;                                   // wave vector [1/m]

    Vector ElementLocalEField(double t, Vector X);
    Vector ElementLocalBField(double t, Vector X);

};

#endif