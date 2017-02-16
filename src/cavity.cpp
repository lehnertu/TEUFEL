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

#include "cavity.h"
#include "global.h"

#include <math.h>

PillboxCavity::PillboxCavity( double freq,
                              double grad,
                              double length,
                              double rcell,
                              int    nlong )
  {
    Freq = freq;
    Grad = grad;
    Length = length;
    R = rcell;
    R2 = R*R;
    Nlong = nlong;
    kz = Nlong*Pi/Length;
    omega = 2.0*Pi*freq;
  };

Vector PillboxCavity::ElementLocalEField(double t, Vector X)
  {
    Vector E = Vector(0.0, 0.0, 0.0);
    if ((X.z>=0.0) && (X.z<=Length))
      {
        double r2 = X.x*X.x + X.y*X.y;
        if (r2<R2)
          {
            // double r = sqrt(r2);
            double Ez = (1-r2/R2)*cos(kz*X.z);
            double Er = 0.5*(1-r2/(2.0*R2))*kz*sin(kz*X.z);     // actually Er/r
            E = Vector(X.x*Er, X.y*Er, Ez) * Grad * cos(omega*t);
          }
      };
    return(E);
  }

Vector PillboxCavity::ElementLocalBField(double t, Vector X)
  {
    return Vector(0.0 ,0.0 ,0.0);
  };
