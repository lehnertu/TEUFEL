/*=========================================================================

  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers

  Module:    externalfield

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

#ifndef EXTERNALFIELD_H
#define EXTERNALFIELD_H

#include <vector>

#include "vector.h"

/*!
    \class ExternalField
    \brief An external field through which to track particles.

    This is the base class for all beam dynamics elements
    treated as external fields (not influenced by the particles)
    during the tracking.

    Methods of this class are called when tracking
    and should be optimized for speed.

    This class takes care of the position and orientation of the elements.
    All methods that are necessary for tracking particles through
    an element are defined here and must be overloaded by
    the classes derived for specific elements.
*/

class ExternalField
{

  public:

    ExternalField();

    virtual ~ExternalField();

    // These functions return the field values
    // queried in global (laboratory frame) coordinates.
    // They take care of the coordinate tranformation
    // and call the appropriate element-specific functions.

    // electric field [MV/m]
    Vector EField(double t, Vector X);

    // magnetic field [T]
    Vector BField(double t, Vector X);

  private:

    // These functions return the activity status and field values
    // queried in element-local coordinates.
    // They are pure virtual function that must be implemented
    // by the particular elements derived from this class.
    virtual Vector ElementLocalEField(double t, Vector X) = 0;
    virtual Vector ElementLocalBField(double t, Vector X) = 0;

};

// ********** Lattice **********

// This is a container class for many objects derived from ExternalField

class Lattice
{

  public:

    // create an empty lattice
    Lattice();

    // the destructor also destroys all
    // elemental external fields belonging to the lattice
    ~Lattice();

    // add some external field to the lattice
    void addElement(ExternalField *element);

    // total electric field [MV/m] of all lattice elements combined
    Vector EField(double t, Vector X);

    // total magnetic field [T] of all lattice elements combined
    Vector BField(double t, Vector X);

private:

    // define an array of pointers to the individual lattice elements
    std::vector<ExternalField *> elements;

};

#endif