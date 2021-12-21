/*=========================================================================

  Program:   TEUFEL - THz Emission from Undulators and Free-Electron Lasers

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
 * Memory is just 6 Bytes for the two Vectors.
 */
class ElMagField
{

public:

    /*! Default constructor:<br>
     * initializes both components to zero
     */
    ElMagField();

    /*! Initializing constructor: */
    ElMagField(Vector E, Vector B);

    /*! Set all components to zero */
    void Zero();

    /*! Electric field report */
    Vector E();

    /*! Magnetic field report */
    Vector B();

    /*! Check if it is different from zero */
    bool isNull() { return (vecE.isNull() and vecB.isNull()); };

    /* Poynting vector - energy flow density */
    Vector Poynting();
    
    /*! Sum of two fields */
    ElMagField operator+ (ElMagField other);

    /*! Accumulating sum of two fields */
    ElMagField& operator+= (ElMagField other);

    /*! Difference of two fields */
    ElMagField operator- (ElMagField other);

    /*! In-place difference of two fields */
    ElMagField& operator-= (ElMagField other);

    /*! Multiplication of the field with a real factor */
    ElMagField operator* (double factor);

    /*! Division of the field by a real factor */
    ElMagField operator/ (double factor);

    /*! In-place multiplication of the field with a real factor */
    ElMagField& operator*= (double factor);

    /*! in-place division of the field by a real factor */
    ElMagField& operator/= (double factor);

private:

    Vector vecE;
    Vector vecB;

};

//! a zero field constant for convenience
const ElMagField ElMagFieldZero(VectorZero,VectorZero);

/*!
 * \class ElMagObs
 * \brief Type for electromagnetic field observations.
 * @author Ulf Lehnert
 * @date 7.4.2017
 * 
 * This class holds the two Vectors E [V/m] and B [T] of the electromagnetic field
 * and in addition the time [s] at which both were observed. This class was created
 * to allow a compact return value of field observations that include the retarded
 * time at which the field arrives at the observation point.
 */
class ElMagObs
{

public:

    /*! Default constructor:<br>
     * initializes both components to zero
     */
    ElMagObs();

    /*! Initializing constructor: */
    ElMagObs(double time, Vector E, Vector B);

    /*! time report */
    double Time();

    /*! Electric field report */
    Vector E();

    /*! Magnetic field report */
    Vector B();

    /*! Electromagnetic field report */
    ElMagField EB();

private:

    double t;
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

    /*! Default destructor */
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
 * \class HomogeneousField
 * \brief A homogeneous electric/magnetic field extending over all space.
 * @author Ulf Lehnert
 * @author Vipul Joshi
 * @date 7.4.2017
 */
class HomogeneousField : public GeneralField
{
    
public:
    
    /*! Default constructor initalizes the fields to zero. */
    HomogeneousField();
    
    /*! Initializing constructor
     * 
     * @param E the electric field strength in V/m
     * @param B the magnetic field in T
     */
    HomogeneousField(Vector E, Vector B);
    
    //! Field report routine
    ElMagField Field(double t, Vector X);
    
private:
    
    ElMagField EB;			// the constant field value
    
};

/*!
 * \class LocalizedField
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
 * TODO: at least we also need a time shift
 */
class LocalizedField : public GeneralField
{

public:

    /*! The default constructor just defines an identity transformation
	of the element coordinates into lab coordinates
    */
    LocalizedField();

    /*! This constructor defines a different origin
     *	of the element local coordinate system, thus, placing the element
     *  at a certain position in lab space and time
     */
    LocalizedField(double time, Vector pos);
    
    /*! All derived classe must provide a destructor */
    virtual ~LocalizedField() {};
    
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

protected :

    /*! This is the position of the element in laboratory frame */
    Vector origin;
    
    /*! This is the local time zero given in lab time */
    double t0;

};

/*!
 * \class Lattice
 * \brief A container class for many objects derived from GeneralField
 * @author Ulf Lehnert
 * @date 7.4.2017
 * 
 * This class holds a number of elements derived from GeneralField.
 * It combines all the fields of the objects it containes to provide particles
 * with the field information needed for tracking. This way the lattice itself
 * is an GeneralField.
 * 
 * Typically the lattice only containes external fields derived from
 * LocalizedField. However, particle bunches may also add the mutual
 * interaction fields of their particles to the lattice.
 * 
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
