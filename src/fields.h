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

#pragma once

#include <vector>
#include "vector.h"

/*!
 * \class ElMagField
 * \brief Type for electromagnetic fields.
 * @author Ulf Lehnert
 * @date 7.4.2017
 * 
 * This class holds the two Vectors E [V/m] and B [T] of the electromagnetic field.
 */
class ElMagField
{

public:
    
    /*! Default constructor:<br>
     * initialized both components to zero
     */
    ElMagField();
    
    /*! Initializing constructor: */
    ElMagField(Vector E, Vector B);
    
    /*! Electric field report */
    Vector E();
    
    /*! Magnetic field report */
    Vector B();

    /*! Sum of two fields */
    ElMagField operator+ (ElMagField other);
    
    /*! Accumulating sum of two fields */
    ElMagField& operator+= (ElMagField other);

private:

    Vector vecE;
    Vector vecB;

};

/*!
 * \class GeneralField
 * \brief An general field through which particles can be tracked.
 * @author Ulf Lehnert
 * @author Vipul Joshi
 * @date 7.4.2017
 * 
 * This class is intended to be a template for all electromagnetic fields
 * handled in this software. Any field derived from this class can be used
 * for tracking particles. It is purely virtual - no implementation is provided.
 */
class GeneralField
{

public:

    /*! Default constructor with no functionality. */
    GeneralField() {};

    virtual ~GeneralField() {};
    /*!
     * The electromagnetic field at a given time and point in space.
     * 
     * The coordinates [m] and the time [s] refer to the laboratory (rest) frame.
     * The field is returned as a tuple of electric field [V/m] and
     * magnetic field [T] vectors.
     * 
     * This method must be overridden by derived classes.
     */
    virtual ElMagField Field(double t, Vector X) = 0;

};

/*!
 * \class ExternalField
 * \brief An external field (lattice element) through which to track particles.
 * @author Ulf Lehnert
 * @author Vipul Joshi
 * @date 7.4.2017
 * 
 * This is the base class for all beam dynamics elements
 * treated as external fields (not influenced by the particles)
 * during the tracking.
 * 
 * This class takes care of the position and orientation of the lattice elements.
 * All elements compute the field in their own (element local) coordinate system.
 * The transformation necessary to compute the field in laboratory coordinates,
 * are handled by this class.
 * 
 * At present only a simple shift of the origin is implemented.
 */
class ExternalField : public GeneralField
{

public:

    /*! The default constructor just defines an identity transformation
	of the element coordinates into lab coordinates
    */
    ExternalField();

    /*! All derived classe must provide a destructor */
    virtual ~ExternalField() {};
    
    /*! The electromagnetic field at a given time and point in space.
     * 
     * The coordinates [m] and the time [s] refer to the laboratory (rest) frame.
     * The field is returned as a tuple of electric field [V/m] and
     * magnetic field [T] vectors.
     * 
     * This calls LocalField() with an appropriately transformed pposition
     * and transforms the returned fields into the laboratory frame.
     */
    virtual ElMagField Field(double t, Vector X);

private:

    /*! The electromagnetic field at a given time and point in space.
     * in element-local coordinates.
     * The coordinates and the time refer to the laboratory (rest) frame.
     * The field is returned as a tuple of electric field [V/m] and
     * magnetic field [T] vectors.
     * 
     * This must be implemented for all derived lattice elements.
     */
    virtual ElMagField LocalField(double t, Vector X) = 0;

private :

    /*! This is the position of the element in laboratory frame */
    Vector origin;

};

/*!
 * \class Lattice
 * \brief This is a container class for many objects derived from GeneralalField
 * @author Ulf Lehnert
 * @date 7.4.2017
 * 
 * This class holds a number of elements derived from GeneralField.
 * It combines all the fields of the objects it containes to provide particles
 * with the field information needed for tracking. This way the lattice itself
 * is an GeneralField.
 * 
 * Typically the lattice only containes external fields derived from
 * ExternalField. However, particle bunches may also add the mutual
 * interaction fields of their particles to the lattice.
 * 
 * \todo The destructor is not yet coded.
 */
class Lattice : public GeneralField
{

public:

    /*! Create an empty lattice. */
    Lattice();

    /*! The destructor also destroys all elements belonging to the lattice */
    ~Lattice();

    /* Add some external field to the lattice */
    void addElement(GeneralField* element);

    /*! Report the number of elements in a lattice */
    int count();

    /*! The electromagnetic field at a given time and point in space
     * of all contained elements combined.
     * The coordinates and the time refer to the laboratory (rest) frame.
     * The field is returned as a tuple of electric field [V/m] and
     * magnetic field [T] vectors.
     */
    ElMagField Field(double t, Vector X);

private:
    
    /*! An array of pointers to the individual lattice elements */
    std::vector<GeneralField*> elements;
    
};
