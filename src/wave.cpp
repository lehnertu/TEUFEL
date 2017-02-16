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

#include "wave.h"
#include "global.h"

#include <math.h>

PlaneWave::PlaneWave( double lambda,
                      double intensity )
  {
    Lambda = lambda;
    omega = 2.0*Pi*SpeedOfLight/Lambda;
    Intensity = 1.0e4 * intensity;      // scale from W/cm² to W/m²
    Epeak = sqrt( 2.0*Intensity / (SpeedOfLight*EpsNull) );
    Bpeak = Epeak / SpeedOfLight;
    a0 = Epeak * SpeedOfLight / mecsquared / omega;
    k = Vector(0.0 ,0.0 ,omega/SpeedOfLight);
  };

double PlaneWave::Omega()
{
  return omega;
}

double PlaneWave::EPeak()
{
  return Epeak;
}

double PlaneWave::BPeak()
{
  return Bpeak;
}

double PlaneWave::A0()
{
  return a0;
}

Vector PlaneWave::ElementLocalEField(double t, Vector X)
  {
    Vector E = Vector(Epeak*cos(omega*t-dot(k,X)), 0.0, 0.0);
    return(E);
  }

Vector PlaneWave::ElementLocalBField(double t, Vector X)
  {
    Vector B = Vector(0.0, Bpeak*cos(omega*t-dot(k,X)), 0.0);
    return(B);
  };
