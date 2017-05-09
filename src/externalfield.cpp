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

#include "externalfield.h"

ExternalField::ExternalField()
{
}

ExternalField::~ExternalField()
{
}

Vector ExternalField::EField(double t, Vector X)
{
  return ElementLocalEField(t, X);
}

Vector ExternalField::BField(double t, Vector X)
{
  return ElementLocalBField(t, X);
}


// ********** Lattice **********

Lattice::Lattice()
{
}

Lattice::~Lattice()
{
}

void Lattice::addElement(ExternalField *element)
{
  // elements->append(element);
  elements.push_back(element);
}

pair<Vector,Vector> Lattice::Field(double t, Vector X)
{
  Vector E = Vector(0.0, 0.0, 0.0);
  Vector B = Vector(0.0, 0.0, 0.0);
  for (unsigned int i=0; i<elements.size(); i++)
  {
    E += elements.at(i)->EField(t, X);
    B += elements.at(i)->BField(t, X);
  };
  return(make_pair(E,B));
}


