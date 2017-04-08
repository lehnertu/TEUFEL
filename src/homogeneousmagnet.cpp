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

#include "homogeneousmagnet.h"
#include "global.h"
#include <math.h>

HomogeneousMagnet::HomogeneousMagnet()
{
    B0=Vector(0.0,0.0,0.0);
}

HomogeneousMagnet::HomogeneousMagnet(Vector B)
{
  B0=B;
}

Vector HomogeneousMagnet::getB0()
{
  return B0;
}

ElMagField HomogeneousMagnet::Field(double t, Vector X)
{
  Vector E=Vector(0.0,0.0,0.0);
  return ElMagField(E,B0);
}
