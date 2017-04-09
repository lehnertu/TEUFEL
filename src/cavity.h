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

#ifndef CAVITY_H
#define CAVITY_H

#include "vector.h"
#include "externalfield.h"

using namespace std;

// a class for a pillbox accelerator cavity

// the cavity is aligned with the z axis
// reference position is the cavity entrance wall
// the mode number is the number of half-wave periods
// in longitudinal direction present inside the cavity - TM(01n)

// the beam direction is along the z axis starting at z=0

/*!
 * \class PillboxCavity
 * \brief Pillbox accelerator cavity.
 * @author Ulf Lehnert
 * @date 7.4.2017
 *
 * \todo update to current field definitions 
 */
class PillboxCavity : public ExternalField
{

  public:

    // constructor
    PillboxCavity( double freq,                        // cavity frequency [Hz]
                   double grad,                        // peak electric field gradient on axis [V/m]
                   double length,                      // length of the cavity [m]
                   double rcell,                       // radius of the cavity [m]
                   int    nlong );                     // mode number in longitudinal direction

  private:

    double Freq;                                // cavity frequency [Hz]
    double Grad;                                // peak electric field gradient on axis [V/m]
    double Length;                              // length of a (half) cell [m]
    double R;                                   // radius of the cavity [m]
    double R2;                                  // the square of it
    int    Nlong;                               // number of (half) cells

    double kz;                                  // cavity periodicity [1/m]
    double omega;                               // cavity frequency [Hz]

    Vector ElementLocalEField(double t, Vector X);
    Vector ElementLocalBField(double t, Vector X);

};

#endif